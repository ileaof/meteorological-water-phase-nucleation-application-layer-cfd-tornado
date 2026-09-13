"""Score the tornado precursor on a sounding, on a saved state, or on both.

    python scripts/tornado_precursor_report.py --sounding
    python scripts/tornado_precursor_report.py --sequence outputs/.../sequence.h5 --time 2790
    python scripts/tornado_precursor_report.py --sounding --sequence ... --time 2790 --csv out.csv

The two halves answer different questions and are printed separately on purpose: the
sounding measures whether the environment can supply the rotation, the state measures
whether the updraft is oriented to tilt it.  See storm_dynamics/tornado_precursor.
"""
from pathlib import Path
import argparse
import csv
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from meteorological_flow.grid import Grid                               # noqa: E402
from storm_dynamics.config import HodographConfig                       # noqa: E402
from storm_dynamics.soundings import build_sounding                     # noqa: E402
from storm_dynamics.tornado_precursor import (environment_precursor,    # noqa: E402
                                              state_precursor)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sounding', action='store_true',
                    help="score this project's analytic supercell sounding")
    ap.add_argument('--sequence', help='a diagnostic sequence.h5 to score a state from')
    ap.add_argument('--time', type=float, help='time in seconds; nearest snapshot is used')
    ap.add_argument('--z-lo', type=float, default=0.0)
    ap.add_argument('--z-hi', type=float, default=1000.0)
    ap.add_argument('--w-updraft', type=float, default=1.0)
    ap.add_argument('--storm-motion', nargs=2, type=float, default=None,
                    metavar=('CX', 'CY'),
                    help='storm motion; omit when the saved fields are already '
                         'storm-relative, as this project\'s sequences are')
    ap.add_argument('--csv', help='append the scored values to this CSV')
    args = ap.parse_args()
    if not args.sounding and not args.sequence:
        ap.error('nothing to score: pass --sounding, --sequence, or both')

    rows = []
    if args.sounding:
        grid = Grid(nx=8, ny=8, nz=48, Lx=72000., Ly=72000., Lz=15000.,
                    z_stretch=1.05, periodic=True)
        base = build_sounding(grid, HodographConfig(kind='quarter_circle', U_max=30.,
                                                    z_turn=3000., u_half=3000.))
        score = environment_precursor(base)
        print(score.summary())
        print()
        rows.append(dict(kind='environment', source='analytic quarter-circle sounding',
                         **score.values,
                         gates_passed=';'.join(score.passed),
                         gates_failed=';'.join(score.failed)))

    if args.sequence:
        import h5py
        with h5py.File(args.sequence, 'r') as f:
            groups = list(f['snapshots'].values())
            times = np.array([g.attrs['time_s'] for g in groups])
            i = 0 if args.time is None else int(np.argmin(np.abs(times - args.time)))
            g = groups[i]
            u, v, w = g['u'][...], g['v'][...], g['w'][...]
            nx, ny, nz = u.shape[0] - 1, v.shape[1] - 1, w.shape[2] - 1
            xc, yc, zf = f['grid/xc'][:], f['grid/yc'][:], f['grid/zf'][:]
        uc = 0.5 * (u[:-1] + u[1:])
        vc = 0.5 * (v[:, :-1] + v[:, 1:])
        wc = 0.5 * (w[:, :, :-1] + w[:, :, 1:])
        grid = Grid(nx=nx, ny=ny, nz=nz, Lx=float(zf[0] + xc[-1] + (xc[1] - xc[0]) / 2),
                    Ly=float(yc[-1] + (yc[1] - yc[0]) / 2), Lz=float(zf[-1]),
                    z_stretch=1.05, periodic=True)
        motion = tuple(args.storm_motion) if args.storm_motion else (0.0, 0.0)
        score = state_precursor(uc, vc, wc, grid, storm_motion=motion,
                                z_lo=args.z_lo, z_hi=args.z_hi, w_updraft=args.w_updraft)
        print('state at t = %.1f s, %s, layer %g-%g m'
              % (times[i], Path(args.sequence).parent.name, args.z_lo, args.z_hi))
        print(score.summary())
        rows.append(dict(kind='state', source=args.sequence, time_s=float(times[i]),
                         z_lo_m=args.z_lo, z_hi_m=args.z_hi, **score.values,
                         gates_passed=';'.join(score.passed),
                         gates_failed=';'.join(score.failed)))

    if args.csv and rows:
        path = Path(args.csv)
        fields = list(dict.fromkeys(k for r in rows for k in r))
        new = not path.exists()
        with path.open('a', newline='', encoding='utf-8') as fh:
            w_ = csv.DictWriter(fh, fieldnames=fields)
            if new:
                w_.writeheader()
            w_.writerows(rows)
        print('\nappended %d row(s) to %s' % (len(rows), path))


if __name__ == '__main__':
    main()
