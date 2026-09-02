"""Command-line interface (ROADMAP §3a):

    python -m atmospheric_data case-info      config/moore_2013.yaml
    python -m atmospheric_data download       config/moore_2013.yaml [--offline]
    python -m atmospheric_data preprocess     config/moore_2013.yaml [--offline] [--outdir DIR]
    python -m atmospheric_data validate-input config/moore_2013.yaml
    python -m atmospheric_data run-case       config/moore_2013.yaml [--steps N]
    python -m atmospheric_data compare-radar  config/moore_2013.yaml

``--offline`` never touches the network.  Exit code is non-zero when QC fails (validate-input)
so the command is CI-usable.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

from .config import CaseConfig
from .cache import Cache


def _load(args):
    cfg = CaseConfig.from_yaml(args.config)
    cfg.offline = bool(getattr(args, "offline", False))
    cache = Cache(cfg.data.cache_directory, offline=cfg.offline)
    return cfg, cache


def _default_outdir(cfg, args):
    return getattr(args, "outdir", None) or os.path.join("outputs", "real_case", cfg.case.name)


def cmd_case_info(args):
    from .sources import synthetic  # noqa: F401 (import check)
    from .sources.base import available
    cfg, _ = _load(args)
    print("=" * 64); print("CASE:", cfg.case.name, "|", cfg.case.date, cfg.case.start_time_utc, "UTC")
    print("=" * 64)
    print(json.dumps(cfg.to_dict(), indent=2))
    print("\nsource availability on this machine:")
    for name, mods in (("HRRR (cfgrib)", ["cfgrib"]), ("ERA5 read (xarray)", ["xarray"]),
                       ("ERA5 download (cdsapi)", ["cdsapi"]), ("NEXRAD (pyart)", ["pyart"]),
                       ("projection (pyproj)", ["pyproj"]), ("sounding/METAR (pandas)", ["pandas"])):
        print("  %-26s %s" % (name, "available" if available(mods) else "MISSING (optional)"))
    return 0


def cmd_download(args):
    cfg, cache = _load(args)
    from .sources import load_atmosphere, load_radar
    try:
        st = load_atmosphere(cfg, cache)
        print("[download] atmosphere ready:", st.ds.attrs.get("source"))
    except Exception as e:
        print("[download] atmosphere:", e)
    if cfg.validation.radar:
        try:
            load_radar(cfg, cache)
        except Exception as e:
            print("[download] radar:", e)
    return 0


def cmd_preprocess(args):
    from . import driver
    cfg, cache = _load(args)
    out = _default_outdir(cfg, args)
    pre = driver.preprocess(cfg, cache, out, max_n=getattr(args, "max_n", 64))
    print("[preprocess] wrote:", out)
    print("  QC report:", pre["qc_md"], "(ok=%s)" % pre["qc"]["summary"]["ok"])
    return 0


def cmd_validate_input(args):
    from . import driver
    cfg, cache = _load(args)
    out = _default_outdir(cfg, args)
    pre = driver.preprocess(cfg, cache, out, max_n=getattr(args, "max_n", 64))
    ok = pre["qc"]["summary"]["ok"]
    print("[validate-input] QC %s (%d/%d)" % ("PASS" if ok else "FAIL",
          pre["qc"]["summary"]["passed"], pre["qc"]["summary"]["total"]))
    return 0 if ok else 2


def cmd_run_case(args):
    from . import driver
    cfg, cache = _load(args)
    out = _default_outdir(cfg, args)
    pre = driver.preprocess(cfg, cache, out, max_n=getattr(args, "max_n", 64))
    if getattr(args, "multilevel", False):
        sims, rep = driver.run_multilevel_real_case(cfg, pre)
        print("[run-case] multilevel cascade: finest dx=%.0f m, zeta=%.3e"
              % (sims[-1].grid.dx, rep["rotation"]["zeta_abs_max"]))
        return 0
    sim = driver.run_case(cfg, pre, steps=getattr(args, "steps", None))
    print("[run-case] done; backend=%s, steps=%d, t=%.1f s" % (sim.grid.backend.name, sim.step, sim.t))
    return 0


def cmd_compare_radar(args):
    from . import driver
    cfg, cache = _load(args)
    out = _default_outdir(cfg, args)
    pre = driver.preprocess(cfg, cache, out, max_n=getattr(args, "max_n", 64))
    sim = driver.run_case(cfg, pre, steps=getattr(args, "steps", None))
    res = driver.compare_radar(cfg, pre, cache, sim=sim)
    with open(os.path.join(out, "radar_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(res["metrics"], f, indent=2)
    print("[compare-radar] metrics ->", os.path.join(out, "radar_metrics.json"))
    return 0


def cmd_storm_watch(args):
    """Operational auto mode: start|status|stop|alerts|cases|retry|replay."""
    from .storm_watch import StormWatchConfig, StormWatchMonitor
    from .storm_watch.db import WatchDB
    action = args.action
    rest = list(args.rest or [])
    # interpret positionals per action: replay/retry take a target first, then optional config
    target, config = None, "config/storm_watch.yaml"
    if action in ("replay", "retry", "map"):
        target = rest[0] if rest else None
        config = rest[1] if len(rest) > 1 else config
    else:
        config = rest[0] if rest else config
    args.target = target
    sw = StormWatchConfig.from_yaml(config) if os.path.exists(config) else StormWatchConfig()
    try:
        base = CaseConfig.from_yaml(config)
    except Exception:
        base = CaseConfig()
    stop_flag = os.path.join(sw.workdir, "STOP")

    if action == "start":
        os.path.exists(stop_flag) and os.remove(stop_flag)
        mon = StormWatchMonitor(sw, base_case_config=base, offline=args.offline, max_n=args.max_n)
        print("[storm-watch] starting (auto_simulate=%s, offline=%s). Ctrl-C or `stop` to end."
              % (sw.actions.auto_simulate, args.offline))
        it = 0
        while args.max_iterations is None or it < args.max_iterations:
            if os.path.exists(stop_flag):
                print("[storm-watch] stop flag seen; exiting"); break
            mon.poll_once(); it += 1
            if args.max_iterations is not None and it >= args.max_iterations:
                break
            import time; time.sleep(sw.alert_poll_seconds)
        mon.close(); return 0
    if action == "stop":
        os.makedirs(sw.workdir, exist_ok=True); open(stop_flag, "w").close()
        print("[storm-watch] stop flag written:", stop_flag); return 0

    db = WatchDB(os.path.join(sw.workdir, "storm_watch.sqlite"))
    if action == "status":
        print(json.dumps({"workdir": sw.workdir, "active_cases": db.active_case_count(),
                          "queued": db.queue_size(), "auto": sw.actions.__dict__}, indent=2))
    elif action == "alerts":
        for a in db.list_alerts():
            print("  %-10s %-28s sev=%-8s status=%s" % (a["level"], a["event"], a["severity"], a["status"]))
    elif action == "cases":
        for c in db.list_cases():
            print("  %s  %-22s  %-22s  (%.2f,%.2f)" % (c["case_id"], c["name"], c["state"],
                                                       c["center_lat"], c["center_lon"]))
    elif action == "replay":
        mon = StormWatchMonitor(sw, base_case_config=base, offline=True, max_n=args.max_n)
        print("[storm-watch] replay:", json.dumps(mon.replay(args.target), indent=2)); mon.close()
    elif action == "retry":
        c = db.get_case(args.target)
        print("[storm-watch] retry", args.target, "->", "not found" if not c else "re-enqueued")
        if c:
            db.set_state(args.target, "DETECTED", "manual retry")
    elif action == "map":
        from .storm_watch import viz
        if not args.target:
            print("[storm-watch] map needs an alert file: storm-watch map ALERT.json [config.yaml]")
            db.close(); return 2
        out = viz.map_alert_file(args.target, cfg=sw)
        print("[storm-watch] map ->", out)
    db.close(); return 0


def main(argv=None):
    p = argparse.ArgumentParser("atmospheric_data", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    swp = sub.add_parser("storm-watch", help="operational auto mode")
    swp.add_argument("action", choices=["start", "status", "stop", "alerts", "cases", "retry", "replay", "map"])
    swp.add_argument("rest", nargs="*", help="[FILE|CASE_ID] and/or config.yaml (see docs)")
    swp.add_argument("--offline", action="store_true")
    swp.add_argument("--max-iterations", type=int, default=None, dest="max_iterations")
    swp.add_argument("--max-n", type=int, default=24, dest="max_n")
    swp.set_defaults(func=cmd_storm_watch)
    cmds = {"case-info": cmd_case_info, "download": cmd_download, "preprocess": cmd_preprocess,
            "validate-input": cmd_validate_input, "run-case": cmd_run_case,
            "compare-radar": cmd_compare_radar}
    for name in cmds:
        sp = sub.add_parser(name)
        sp.add_argument("config", help="path to the case YAML")
        sp.add_argument("--offline", action="store_true", help="never access the network")
        sp.add_argument("--outdir", default=None, help="output directory")
        sp.add_argument("--max-n", type=int, default=64, dest="max_n",
                        help="cap grid points per axis (test/dev; raise for production)")
        sp.add_argument("--steps", type=int, default=None, help="run-case: number of steps")
        sp.add_argument("--multilevel", action="store_true",
                        help="run-case: drive the AMR parent->nest->fine cascade from the real IC")
    args = p.parse_args(argv)
    if hasattr(args, "func"):
        return args.func(args)
    return cmds[args.command](args)
