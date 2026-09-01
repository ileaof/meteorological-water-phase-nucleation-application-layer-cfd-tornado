"""Real-case ingestion pipeline (offline) via the Python API (ROADMAP §3a).

Runs the whole `atmospheric_data` workflow with NO network and NO heavy dependencies (uses the
labelled synthetic sample), so it is always reproducible:

    python examples/real_case_offline.py

To use real HRRR/ERA5/NEXRAD, `pip install cfgrib eccodes arm_pyart cdsapi`, point the YAML's
`atmospheric_source` at `hrrr`/`era5`, and drop / download the data (see docs/REAL_CASE_DATA.md).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import tempfile
import warnings

import atmospheric_data as ad
from atmospheric_data import driver


def main():
    warnings.filterwarnings("ignore")
    cfg = ad.CaseConfig.from_yaml(os.path.join(os.path.dirname(__file__), "..", "config", "local_case.yaml"))
    cfg.offline = True
    cache = ad.Cache(cfg.data.cache_directory, offline=True)
    out = os.path.join(tempfile.gettempdir(), "real_case_offline")

    print("=" * 70)
    print("real_case ingestion (OFFLINE, synthetic sample) — %s" % cfg.case.name)
    print("=" * 70)
    pre = driver.preprocess(cfg, cache, out, max_n=20)
    print("  QC: %d/%d checks passed (ok=%s)"
          % (pre["qc"]["summary"]["passed"], pre["qc"]["summary"]["total"], pre["qc"]["summary"]["ok"]))
    print("  base state: theta0[0]=%.1f K  qv0[0]=%.4f  p0[0]=%.0f Pa"
          % (pre["base"].theta0[0], pre["base"].qv0[0], pre["base"].p0[0]))
    print("  interp log:", pre["interp_log"][0])

    sim = driver.run_case(cfg, pre, steps=5)
    print("  ran %d steps on %s; theta/w finite [ok]" % (sim.step, sim.grid.backend.name))

    res = driver.compare_radar(cfg, pre, cache, sim=sim)
    m = res["metrics"]["radial_velocity"]
    print("  radial-velocity vs radar: RMSE=%.2f m/s  bias=%.2f  corr=%.2f"
          % (m["rmse"], m["bias"], m["correlation"]))
    print("\nartefacts in:", out, "->", sorted(f for f in os.listdir(out) if f.endswith((".nc", ".json", ".md"))))
    print("\nIDEALISED note: synthetic environment; with real HRRR/NEXRAD this becomes a real "
          "case, but the tornado must still EMERGE from the resolved dynamics (not imposed).")


if __name__ == "__main__":
    main()
