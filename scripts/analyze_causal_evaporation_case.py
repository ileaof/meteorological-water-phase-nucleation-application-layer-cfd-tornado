"""Generate the invariant single-case diagnostics used by the evaporation pilot.

This deliberately omits the historical-sequence prose report: the causal pilot
is a restart experiment and its interpretation belongs in the cross-case report.
"""
from pathlib import Path
import argparse
import hashlib
import json

import h5py
import numpy as np

from analyze_diagnostic_sequence import (
    ROOT,
    audit_and_track,
    lagrangian_budget,
    trajectories,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sequence", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seed-time", type=float, default=2400.0)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    source_files = [
        ROOT / "scripts" / "analyze_diagnostic_sequence.py",
        Path(__file__).resolve(),
        ROOT / "src" / "storm_dynamics" / "diagnostic_capture.py",
    ]
    with h5py.File(args.sequence, "r") as f:
        summary, peak_index, seeds, track = audit_and_track(
            f, args.out, seed_time=args.seed_time
        )
        positions = trajectories(f, peak_index, seeds, max_step=2.0)
        fine = trajectories(f, peak_index, seeds, max_step=1.0)
        coarse = trajectories(f, peak_index, seeds, max_step=2.0, stride=2)
        distances = lambda other: np.concatenate(
            [np.linalg.norm(positions[i] - p, axis=1) for i, p in other.items()]
        )
        d_fine = distances(fine)
        d_coarse = distances(coarse)
        sensitivity = dict(
            rk4_2s_vs_1s_max_m=float(np.nanmax(d_fine)),
            rk4_2s_vs_1s_median_m=float(np.nanmedian(d_fine)),
            cadence_30_vs_60s_max_m=float(np.nanmax(d_coarse)),
            cadence_30_vs_60s_median_m=float(np.nanmedian(d_coarse)),
            finite_pairs_rk4=int(np.isfinite(d_fine).sum()),
            finite_pairs_cadence=int(np.isfinite(d_coarse).sum()),
        )
        result = lagrangian_budget(f, positions, args.out)
        summary.update(
            input_sequence=str(args.sequence.resolve()),
            input_sequence_sha256=hashlib.sha256(args.sequence.read_bytes()).hexdigest(),
            analysis_source_sha256={
                str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in source_files
            },
            trajectory_sensitivity=sensitivity,
            lagrangian={
                key: value for key, value in result.items()
                if key not in ("records", "budget")
            },
            restart_experiment=True,
            intervention_start_s=float(f.attrs.get("restart_time_s", 2370.211552)),
            stored_vertical_extent_m=float(f["grid/zc"][int(f.attrs["budget_nz_with_halo"]) - 1]),
        )
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8"
    )
    print(f"Complete: {args.out}", flush=True)


if __name__ == "__main__":
    main()
