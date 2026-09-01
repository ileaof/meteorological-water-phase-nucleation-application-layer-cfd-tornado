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


def main(argv=None):
    p = argparse.ArgumentParser("atmospheric_data", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
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
    args = p.parse_args(argv)
    return cmds[args.command](args)
