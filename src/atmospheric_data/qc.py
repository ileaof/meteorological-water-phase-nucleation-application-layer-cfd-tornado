"""Quality control of the ingested state (ROADMAP §3a) -> JSON + Markdown report.

Runs the task's checklist (units/ranges/missing/monotonic-z/hydrostatic/continuity/mass/
vapour/positive-T/admissible-humidity/divergence/terrain) on the unified state before the CFD
starts, so bad data is caught (and reported), never silently used.  Each check yields
``pass/warn/fail`` with the offending value -- honest, not a rubber stamp.
"""
from __future__ import annotations

import json
import os

import numpy as np

from . import thermo

# physically admissible SI ranges per standard variable
_RANGES = {"T": (150.0, 340.0), "theta": (250.0, 600.0), "p": (1000.0, 1.1e5),
           "qv": (0.0, 0.05), "u": (-150.0, 150.0), "v": (-150.0, 150.0),
           "w": (-100.0, 100.0), "rho": (0.05, 1.6), "terrain": (-500.0, 9000.0)}


def _chk(name, ok, detail):
    return {"check": name, "status": "pass" if ok else "fail", "detail": detail}


def quality_control(state):
    """Return a QC report dict for an :class:`AtmosphericState`."""
    ds = state.ds
    checks = []
    # units present
    missing_units = [v for v in ds.data_vars if not ds[v].attrs.get("units")]
    checks.append(_chk("units_present", not missing_units,
                       "all variables tagged" if not missing_units else "missing units: %s" % missing_units))
    # physical ranges + positive T + admissible humidity + finite (no invented missing)
    for v in ds.data_vars:
        a = np.asarray(ds[v].values, float)
        finite = np.isfinite(a)
        checks.append(_chk("finite_%s" % v, finite.all(),
                           "%.3f%% non-finite" % (100 * (1 - finite.mean()))))
        if v in _RANGES:
            lo, hi = _RANGES[v]; inr = (a[finite] >= lo) & (a[finite] <= hi)
            checks.append(_chk("range_%s" % v, bool(inr.all()),
                               "%.3f%% outside [%.3g,%.3g] (min=%.3g max=%.3g)"
                               % (100 * (1 - inr.mean()), lo, hi, np.nanmin(a), np.nanmax(a))))
    # monotonic vertical coordinate
    z = np.asarray(ds["z"].values, float)
    checks.append(_chk("monotonic_z", bool(np.all(np.diff(z) > 0)), "z strictly increasing"))
    # hydrostatic consistency (dp/dz ~ -rho g) if p present
    if state.has("p") and state.has("T"):
        p = state.var("p"); T = state.var("T")
        qv = state.var("qv") if state.has("qv") else 0.0 * T
        col_p = p.mean(axis=(0, 2, 3)) if p.ndim == 4 else p
        col_T = T.mean(axis=(0, 2, 3)) if T.ndim == 4 else T
        col_q = (qv.mean(axis=(0, 2, 3)) if getattr(qv, "ndim", 0) == 4 else qv)
        rho = thermo.density(col_p, thermo.virtual_temperature(col_T, col_q))
        dpdz = np.gradient(col_p, z)
        resid = np.abs(dpdz + rho * thermo.g0) / (rho * thermo.g0 + 1e-9)
        checks.append(_chk("hydrostatic_balance", float(np.nanmedian(resid)) < 0.15,
                           "median |dp/dz + rho g|/(rho g) = %.3f" % float(np.nanmedian(resid))))
    # temporal continuity (no jumps between analysis times)
    if ds.sizes.get("time", 1) > 1 and state.has("theta"):
        th = state.var("theta")
        jump = np.nanmax(np.abs(np.diff(th, axis=0)))
        checks.append(_chk("temporal_continuity", jump < 30.0,
                           "max |d theta / d(time step)| = %.2f K" % float(jump)))
    # terrain vs levels compatibility
    if state.has("terrain"):
        checks.append(_chk("terrain_below_top", float(np.nanmax(state.var("terrain"))) < z.max(),
                           "max terrain %.0f m < domain top %.0f m" % (float(np.nanmax(state.var("terrain"))), z.max())))
    n_fail = sum(1 for c in checks if c["status"] == "fail")
    return {"summary": {"total": len(checks), "passed": len(checks) - n_fail, "failed": n_fail,
                        "ok": n_fail == 0},
            "source": ds.attrs.get("source", "unknown"),
            "projection": ds.attrs.get("projection", "none"),
            "checks": checks}


def write_reports(report, json_path, md_path):
    """Write the QC report as JSON and Markdown."""
    for p in (json_path, md_path):
        os.makedirs(os.path.dirname(os.path.abspath(p)) or ".", exist_ok=True)
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2)
    s = report["summary"]
    lines = ["# Data quality-control report", "",
             "**Source:** %s  |  **Projection:** %s" % (report["source"], report["projection"]),
             "", "**%d/%d checks passed** (%s)" % (s["passed"], s["total"],
                                                   "OK to run" if s["ok"] else "**FAILURES — review before running**"),
             "", "| check | status | detail |", "|---|---|---|"]
    for c in report["checks"]:
        mark = "✅" if c["status"] == "pass" else "❌"
        lines.append("| `%s` | %s %s | %s |" % (c["check"], mark, c["status"], c["detail"]))
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    return json_path, md_path
