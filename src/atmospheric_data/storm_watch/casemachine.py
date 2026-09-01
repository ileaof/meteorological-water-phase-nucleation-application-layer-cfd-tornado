"""Case state machine (storm-watch).

DETECTED → DOWNLOADING → READY_FOR_PREPROCESSING → PREPROCESSING → READY_FOR_SIMULATION
        → (SIMULATING → VALIDATING → COMPLETED)   [only if auto_simulate AND the trigger is met]
        → FAILED (on error; the monitor retries with backoff)

Builds a real_case :class:`atmospheric_data.CaseConfig` from the alert + auto domain, drives the
existing ``driver`` pipeline, and logs every transition to SQLite. Auto-simulation is **off by
default** and gated by ``simulation_trigger`` (min level/severity).
"""
from __future__ import annotations

import copy
import hashlib
import os

from ..config import CaseConfig
from .. import driver
from ..cache import Cache
from . import domain as dom

_LEVEL_RANK = {"watch": 1, "warning": 2, "confirmed": 3}
_SEV_RANK = {"unknown": 0, "minor": 1, "moderate": 2, "severe": 3, "extreme": 4}
STATES = ("DETECTED", "DOWNLOADING", "READY_FOR_PREPROCESSING", "PREPROCESSING",
          "READY_FOR_SIMULATION", "SIMULATING", "VALIDATING", "COMPLETED", "FAILED")


def case_id_for(alert):
    """Stable, unique case id from the alert id (dedup across restarts)."""
    return "case_" + hashlib.sha1(alert.alert_id.encode("utf-8")).hexdigest()[:12]


def _trigger_ok(level, severity, sw):
    st = sw.simulation_trigger
    return (sw.actions.auto_simulate
            and _LEVEL_RANK.get(level, 0) >= _LEVEL_RANK.get(st.minimum_level, 99)
            and _SEV_RANK.get(severity, 0) >= _SEV_RANK.get(st.minimum_severity, 99))


def build_case_config(alert, level, sw, base=None, offline=False):
    """A real_case CaseConfig for this alert (domain from the polygon, nearest radar)."""
    cfg = copy.deepcopy(base) if base is not None else CaseConfig()
    d = dom.build_domain(alert, sw)
    radars = dom.nearest_radars(d["center_lat"], d["center_lon"], sw)
    cfg.case.name = "%s_%s" % (alert.event.replace(" ", "_").lower(), case_id_for(alert)[5:])
    onset = (alert.onset_time or alert.effective_time or "2013-05-20T20:00:00Z")
    cfg.case.date = onset[:10]
    cfg.case.start_time_utc = onset[11:16] if len(onset) >= 16 else "20:00"
    cfg.domain.center_lat = d["center_lat"]; cfg.domain.center_lon = d["center_lon"]
    cfg.domain.width_km = d["width_km"]; cfg.domain.height_km = d["height_km"]
    cfg.data.atmospheric_source = "hrrr"; cfg.data.fallback_source = "era5"
    cfg.data.radar_station = radars[0][0] if radars else cfg.data.radar_station
    cfg.model.input_mode = "real_case"
    cfg.offline = offline
    cfg.validate()
    return cfg, d, radars


def advance_case(alert, level, sw, db, base=None, offline=False, logger=print, max_n=24, do_simulate=None):
    """Run the case through the state machine.  Returns the final state string.  Stops at
    READY_FOR_SIMULATION unless auto-simulate is enabled and the trigger is met."""
    cid = case_id_for(alert)
    workdir = os.path.join(sw.workdir, cid)
    os.makedirs(workdir, exist_ok=True)
    cfg, dom_info, radars = build_case_config(alert, level, sw, base=base, offline=offline)
    hrrr_run = dom.latest_hrrr_run(alert.onset_time or alert.effective_time or cfg.case.date + "T20:00:00Z")
    metars = dom.nearest_metars(cfg.domain.center_lat, cfg.domain.center_lon)
    db.upsert_case(cid, alert.alert_id, cfg.case.name, "DETECTED", alert.severity, level,
                   cfg.domain.center_lat, cfg.domain.center_lon, workdir)
    logger("[case %s] DETECTED %s (level=%s sev=%s) domain=(%.2f,%.2f,%.0fkm) radars=%s HRRR=%s"
           % (cid, alert.event, level, alert.severity, cfg.domain.center_lat, cfg.domain.center_lon,
              cfg.domain.width_km, [r[0] for r in radars], hrrr_run))
    cache = Cache(cfg.data.cache_directory, offline=offline)
    try:
        db.set_state(cid, "DOWNLOADING", "auto_download=%s" % sw.actions.auto_download)
        if sw.actions.auto_download:
            from ..sources import load_atmosphere, load_radar
            load_atmosphere(cfg, cache, logger=lambda *a: None)
            if sw.data.nexrad:
                load_radar(cfg, cache, logger=lambda *a: None)
        db.set_state(cid, "READY_FOR_PREPROCESSING")

        if not sw.actions.auto_preprocess:
            return "READY_FOR_PREPROCESSING"
        db.set_state(cid, "PREPROCESSING")
        pre = driver.preprocess(cfg, cache, workdir, logger=lambda *a: None, max_n=max_n)
        ok = pre["qc"]["summary"]["ok"]
        db.set_state(cid, "READY_FOR_SIMULATION", "QC ok=%s (%d/%d)"
                     % (ok, pre["qc"]["summary"]["passed"], pre["qc"]["summary"]["total"]))
        logger("[case %s] READY_FOR_SIMULATION -- IC/BC + QC (%s), NetCDF in %s"
               % (cid, "ok" if ok else "REVIEW", workdir))

        run_sim = _trigger_ok(level, alert.severity, sw) if do_simulate is None else do_simulate
        if not run_sim:
            return "READY_FOR_SIMULATION"
        db.set_state(cid, "SIMULATING", "trigger met (auto_simulate)")
        sim = driver.run_case(cfg, pre, steps=5, logger=lambda *a: None)
        db.set_state(cid, "VALIDATING")
        driver.compare_radar(cfg, pre, cache, sim=sim, logger=lambda *a: None)
        db.set_state(cid, "COMPLETED")
        logger("[case %s] COMPLETED" % cid)
        return "COMPLETED"
    except Exception as e:                                 # pragma: no cover (error path)
        db.set_state(cid, "FAILED", note=str(e)[:200]); db.bump_retry(cid)
        logger("[case %s] FAILED: %s" % (cid, e))
        return "FAILED"
