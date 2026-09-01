"""storm-watch operational auto mode tests (ROADMAP §3a).

Exercises the full automation offline (no network, no heavy deps) via alert replay and crafted
alerts: detection, dedup, filtering, level classification (warning vs confirmed), auto-domain,
radar/HRRR selection, the case state machine, resource limits, queue persistence, API-failure
resilience, and the safe auto_simulate default.
"""
import os
import tempfile
import warnings

import numpy as np  # noqa: F401
import pytest

from atmospheric_data.storm_watch import (StormWatchConfig, WatchDB, StormWatchMonitor,
                                          alerts as A, domain as D, casemachine)

SAMPLE = os.path.join(os.path.dirname(__file__), "data", "sample_tornado_alert.json")


def _sw(tmp):
    sw = StormWatchConfig()
    sw.workdir = tmp
    sw.notifications.desktop = False
    return sw


def _alert(event="Tornado Warning", sev="extreme", desc="", poly=True, aid="a1"):
    return A.Alert(alert_id=aid, event=event, severity=sev, description=desc,
                   onset_time="2013-05-20T20:00:00Z",
                   polygon=[[-97.6, 35.3], [-97.3, 35.3], [-97.3, 35.5], [-97.6, 35.5]] if poly else [])


# ---- classification: warning vs confirmed (scientific safety) -------------------
def test_classify_warning_vs_confirmed():
    sw = StormWatchConfig()
    plain = _alert(desc="a tornado warning has been issued, take cover")
    assert A.classify_level(plain, sw) == "warning"                 # NOT confirmed
    confirmed = _alert(desc="this is a RADAR CONFIRMED TORNADO with tornado debris signature")
    assert A.classify_level(confirmed, sw) == "confirmed"
    assert A.confirmation_stage(confirmed) in ("radar-indicated", "observed", "damage-confirmed")
    svr = _alert(event="Severe Thunderstorm Warning", sev="severe", desc="60 mph winds")
    assert A.classify_level(svr, sw) == "warning"
    filtered = _alert(event="Flood Advisory", sev="minor")
    assert A.classify_level(filtered, sw) is None                   # event/severity filtered out


# ---- auto domain + asset selection ---------------------------------------------
def test_auto_domain_radars_and_hrrr():
    sw = StormWatchConfig()
    al = _alert(desc="moving northeast at 30 mph")
    dom = D.build_domain(al, sw)
    assert 34 < dom["center_lat"] < 37 and -99 < dom["center_lon"] < -96
    assert dom["width_km"] == max(sw.automatic_domain.upstream_margin_km + sw.automatic_domain.downstream_margin_km,
                                  2 * sw.automatic_domain.lateral_margin_km)
    radars = D.nearest_radars(dom["center_lat"], dom["center_lon"], sw)
    assert radars and radars[0][0] in D.NEXRAD_STATIONS and len(radars) <= sw.radar.maximum_radars
    day, hh = D.latest_hrrr_run("2013-05-20T20:12:00Z")
    assert day == "20130520" and hh == "20"


def test_missing_polygon_raises_in_domain():
    with pytest.raises(ValueError):
        D.build_domain(_alert(poly=False), StormWatchConfig())


# ---- state machine reaches READY_FOR_SIMULATION (auto_simulate off) --------------
def test_case_reaches_ready_for_simulation_offline():
    tmp = tempfile.mkdtemp(); sw = _sw(tmp)
    db = WatchDB(os.path.join(tmp, "w.sqlite"))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        state = casemachine.advance_case(_alert(), "confirmed", sw, db, offline=True,
                                         logger=lambda *a: None, max_n=16)
    assert state == "READY_FOR_SIMULATION"                          # stops here (auto_simulate=False)
    cid = casemachine.case_id_for(_alert())
    tr = [t["to_state"] for t in db.transitions(cid)]
    assert tr[:5] == ["DETECTED", "DOWNLOADING", "READY_FOR_PREPROCESSING", "PREPROCESSING", "READY_FOR_SIMULATION"]
    assert os.path.exists(os.path.join(sw.workdir, cid, "initial_conditions.nc"))  # internal NetCDF
    db.close()


# ---- dedup + monitor restart persistence ----------------------------------------
def test_dedup_and_restart_persistence():
    tmp = tempfile.mkdtemp(); sw = _sw(tmp)
    db_path = os.path.join(tmp, "w.sqlite")
    mon = StormWatchMonitor(sw, db_path=db_path, offline=True, logger=lambda *a: None, max_n=16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s1 = mon.poll_once([_alert(aid="dup1")])
        s2 = mon.poll_once([_alert(aid="dup1")])                    # same alert again
    assert s1["new"] == 1 and s2["new"] == 0 and s2["skipped"] == 1
    mon.close()
    mon2 = StormWatchMonitor(sw, db_path=db_path, offline=True, logger=lambda *a: None)  # "restart"
    assert mon2.db.alert_seen("dup1")                              # persisted across restart
    mon2.close()


# ---- two simultaneous alerts + replay (completion criterion) --------------------
def test_replay_two_alerts_creates_two_cases():
    tmp = tempfile.mkdtemp(); sw = _sw(tmp)
    mon = StormWatchMonitor(sw, db_path=os.path.join(tmp, "w.sqlite"), offline=True,
                            logger=lambda *a: None, max_n=16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        summary = mon.replay(SAMPLE)
    assert summary["new"] == 2 and summary["advanced"] == 2
    cases = mon.db.list_cases()
    assert len(cases) == 2 and all(c["state"] == "READY_FOR_SIMULATION" for c in cases)
    levels = {a["level"] for a in mon.db.list_alerts()}
    assert "confirmed" in levels and "warning" in levels           # tornado=confirmed, svr=warning
    mon.close()


# ---- resource limit: 3rd alert queued, not advanced -----------------------------
def test_resource_limit_queues_extra_cases():
    tmp = tempfile.mkdtemp(); sw = _sw(tmp)
    sw.resource_limits.maximum_active_cases = 1
    mon = StormWatchMonitor(sw, db_path=os.path.join(tmp, "w.sqlite"), offline=True,
                            logger=lambda *a: None, max_n=16)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mon.poll_once([_alert(aid="c1"), _alert(aid="c2"), _alert(aid="c3")])
    assert mon.db.active_case_count() >= 1 and mon.db.queue_size() >= 1   # extras queued
    assert mon.db.dequeue() is not None                              # queue persists / drains
    mon.close()


# ---- API failure resilience -----------------------------------------------------
def test_api_failure_does_not_crash_monitor(monkeypatch):
    tmp = tempfile.mkdtemp(); sw = _sw(tmp)
    mon = StormWatchMonitor(sw, db_path=os.path.join(tmp, "w.sqlite"), offline=False,
                            logger=lambda *a: None)
    def boom(*a, **k):
        raise RuntimeError("api down")
    monkeypatch.setattr("atmospheric_data.storm_watch.alerts.fetch_active_alerts", boom)
    out = mon.poll_once()                                            # fetch fails
    assert out["fetched"] == 0 and out["fail_streak"] == 1          # logged, did not raise
    mon.close()


# ---- auto_simulate default is OFF (safety) --------------------------------------
def test_auto_simulate_default_off():
    assert StormWatchConfig().actions.auto_simulate is False
