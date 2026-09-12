"""Case-control test for the 2370-2610 s "preparatory phase" LES dominance.

The 27 published parcels were seeded AT the 2790 s vortex peak and integrated
backward, so every statement of the form "positive in 27/27 parcels" is
conditioned on the outcome.  This script builds the missing control arm.

PRE-REGISTERED SEEDING RULE (fixed before any result was inspected)
-------------------------------------------------------------------
Seeds are placed at t = 1800.165633 s (snapshot 60), the START of the analysis
window, and integrated FORWARD.  Selection therefore cannot know the outcome.

  levels     zc[1..3] = 121.66, 207.52, 297.67 m  (brackets the treatment's own
             1800 s height range of about 60-310 m)
  near       every 2nd grid column (1200 m) with r <= 9 km of the treatment
             centroid at 1800 s; the treatment's own 1800 s spread is 6.3 km
  far        every 4th grid column (2400 m) with 9 km < r <= 18 km
  excluded   2 columns at each domain border (the tracker's own rule)

No filtering on zeta, omega_h, w, buoyancy or anything else.  Parcels leaving
the 0-2000 m budget volume or the domain are marked invalid per interval,
exactly as in the published analysis.

The OUTCOME (distance to the tracked low-level peak at 2790 s, and zeta there)
is measured only after the budgets are computed, and is never used to select.

ARMS
----
  validate   the treatment's own recorded positions, re-sampled by this code;
             diffed against analysis/lagrangian_budget.csv to prove the two
             budget samplers agree before any control number is compared
  retrace    the treatment's 1800 s positions integrated FORWARD to 2850 s;
             bounds the backward/forward integration asymmetry
  control    the pre-registered lattice above

Reads only.  No simulation is run, no existing output is modified.
"""
from pathlib import Path
import argparse
import csv
import json
import sys

import numpy as np
import h5py

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from storm_dynamics.diagnostic_capture import curl                              # noqa: E402
from analyze_diagnostic_sequence import VelocityInterval, rk4, native, write_csv  # noqa: E402
from scipy.interpolate import RegularGridInterpolator                           # noqa: E402

PHASES = [(60, 79), (79, 87), (87, 93), (93, 95)]
SEED_INDEX = 60
LAST_INDEX = 95
NK = None


def seed_lattice(coords, centre, zlevels, near_stride=2, far_stride=4,
                 near_radius=9000.0, far_radius=18000.0):
    """Deterministic lattice; no property of the flow enters the choice."""
    xc, yc = coords[0], coords[1]
    seeds, strata = [], []
    for i in range(2, len(xc) - 2):
        for j in range(2, len(yc) - 2):
            r = float(np.hypot(xc[i] - centre[0], yc[j] - centre[1]))
            if r <= near_radius and i % near_stride == 0 and j % near_stride == 0:
                stratum = 'near'
            elif near_radius < r <= far_radius and i % far_stride == 0 and j % far_stride == 0:
                stratum = 'far'
            else:
                continue
            for z in zlevels:
                seeds.append((xc[i], yc[j], z))
                strata.append((stratum, r))
    return np.array(seeds, float), strata


def integrate_forward(f, groups, seeds, i0, i1):
    positions = {i0: seeds.copy()}
    current = seeds.copy()
    for i in range(i0, i1):
        interval = VelocityInterval(groups[i], groups[i + 1], f)
        current = rk4(interval, current, interval.t0, interval.t1, 2.0)
        positions[i + 1] = current.copy()
        print('  trajectory %d/%d' % (i + 1, i1), flush=True)
    return positions


def budgets(f, groups, coords, positions, i0, i1):
    """Midpoint-quadrature material budget, identical in form to the published one."""
    def sample(a, p):
        return RegularGridInterpolator(coords, a, bounds_error=False, fill_value=np.nan)(p)

    n = len(positions[i0])
    keys = None
    per_interval = {}
    previous = None
    for i in range(i0, i1 + 1):
        g = groups[i]
        om = curl(native(g, NK), coords)
        values = np.stack([sample(a, positions[i]) for a in om], axis=-1)
        if previous is not None:
            mid = (positions[i] + positions[i - 1]) / 2
            inside = ((mid[:, 2] <= 2000) & (positions[i][:, 2] <= 2000)
                      & (positions[i - 1][:, 2] <= 2000))
            ki = g['kinematic_integrals']
            prod = {}
            for name in ('vector_stretching', 'vector_dilatation'):
                prod[name] = np.stack([sample(a, mid) for a in ki[name][:]], axis=-1)
            axial = sample(ki['stretching'][:], mid)
            tilt = sample(ki['tilting'][:], mid)
            adv = curl(native(g['increments/advection']), coords)
            continuous = sum(ki[name][:] for name in
                             ('vector_advection', 'vector_stretching', 'vector_dilatation'))
            prod['advection_operator_remainder'] = np.stack(
                [sample(a, mid) for a in adv - continuous], axis=-1)
            for name, stage in g['increments'].items():
                if name != 'advection':
                    prod[name] = np.stack([sample(a, mid) for a in curl(native(stage), coords)],
                                          axis=-1)
            total = sum(prod.values())
            change = values - previous
            if keys is None:
                keys = (['delta_zeta', 'residual_zeta', 'stretching_zeta', 'tilting_zeta']
                        + [name + '_zeta' for name in sorted(prod)])
            row = {}
            row['delta_zeta'] = change[:, 2]
            row['residual_zeta'] = change[:, 2] - total[:, 2]
            row['stretching_zeta'] = axial
            row['tilting_zeta'] = tilt
            for name, val in prod.items():
                row[name + '_zeta'] = val[:, 2]
            row['valid'] = inside & np.isfinite(total[:, 2]) & np.isfinite(change[:, 2])
            per_interval[i] = row
        previous = values
        print('  budget %d/%d' % (i, i1), flush=True)
    return per_interval, keys, n


def phase_integrals(per_interval, keys, n):
    out = []
    for a, b in PHASES:
        acc = {k: np.zeros(n) for k in keys}
        valid = np.ones(n, bool)
        for i in range(a + 1, b + 1):
            row = per_interval[i]
            valid &= row['valid']
            for k in keys:
                acc[k] += np.nan_to_num(row[k])
        out.append((a, b, acc, valid))
    return out


def main():
    global NK
    ap = argparse.ArgumentParser()
    ap.add_argument('--sequence', default='outputs/diagnostic_sequence_20260905/sequence.h5')
    ap.add_argument('--analysis', default='outputs/diagnostic_sequence_20260905/analysis')
    ap.add_argument('--out', default='outputs/control_parcels_20260912')
    # The defaults are the pre-registered rule.  Overrides only increase the
    # sample on the SAME rule; they never change what is selected on.
    ap.add_argument('--near-stride', type=int, default=2)
    ap.add_argument('--far-stride', type=int, default=4)
    ap.add_argument('--near-radius', type=float, default=9000.0)
    ap.add_argument('--far-radius', type=float, default=18000.0)
    ap.add_argument('--levels', default='1,2,3',
                    help='zc indices to seed at; default brackets the treatment')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    analysis = Path(args.analysis)

    f = h5py.File(args.sequence, 'r')
    groups = list(f['snapshots'].values())
    times = np.array([g.attrs['time_s'] for g in groups])
    NK = int(f.attrs['budget_nz_with_halo'])
    coords = tuple(f['grid'][n][:] for n in ('xc', 'yc', 'zc'))
    coords = (coords[0], coords[1], coords[2][:NK])

    tre = {}
    for r in csv.DictReader((analysis / 'parcels.csv').open()):
        i = int(r['index'])
        if SEED_INDEX <= i <= LAST_INDEX:
            tre.setdefault(i, {})[int(r['parcel'])] = (float(r['x_m']), float(r['y_m']),
                                                       float(r['z_m']))
    treatment = {i: np.array([tre[i][p] for p in sorted(tre[i])]) for i in tre}
    ntre = len(treatment[SEED_INDEX])
    centre = treatment[SEED_INDEX][:, :2].mean(axis=0)
    print('treatment centroid at 1800 s: %.1f, %.1f m; n=%d' % (centre[0], centre[1], ntre),
          flush=True)

    summary = {
        'seed_time_s': float(times[SEED_INDEX]),
        'last_time_s': float(times[LAST_INDEX]),
        'treatment_centroid_1800_m': [float(centre[0]), float(centre[1])],
        'treatment_n': ntre,
        'phase_boundaries_s': [float(times[a]) for a, _ in PHASES]
                              + [float(times[PHASES[-1][1]])],
    }

    print('ARM retrace: trajectories', flush=True)
    retrace_pos = integrate_forward(f, groups, treatment[SEED_INDEX], SEED_INDEX, LAST_INDEX)
    err = np.linalg.norm(retrace_pos[LAST_INDEX] - treatment[LAST_INDEX], axis=1)
    summary['retrace_endpoint_error_m'] = dict(median=float(np.nanmedian(err)),
                                               max=float(np.nanmax(err)))

    print('ARM control: trajectories', flush=True)
    zlev = [float(coords[2][int(k)]) for k in args.levels.split(',')]
    seeds, strata = seed_lattice(coords, centre, zlev, args.near_stride, args.far_stride,
                                 args.near_radius, args.far_radius)
    summary['seed_rule'] = dict(levels_m=zlev, near_stride=args.near_stride,
                                far_stride=args.far_stride, near_radius_m=args.near_radius,
                                far_radius_m=args.far_radius)
    print('  %d control seeds' % len(seeds), flush=True)
    ctrl_pos = integrate_forward(f, groups, seeds, SEED_INDEX, LAST_INDEX)

    # One budget pass over the snapshots for all three arms at once; the file is
    # 13 GB and each pass re-reads every operator increment.
    nctrl = len(seeds)
    combined = {i: np.vstack([treatment[i], retrace_pos[i], ctrl_pos[i]])
                for i in range(SEED_INDEX, LAST_INDEX + 1)}
    print('budgets for %d parcels (%d treatment + %d retrace + %d control)'
          % (len(combined[SEED_INDEX]), ntre, ntre, nctrl), flush=True)
    allper, keys, _ = budgets(f, groups, coords, combined, SEED_INDEX, LAST_INDEX)

    def slice_arm(lo, hi):
        return {i: {k: (v[lo:hi] if hasattr(v, '__len__') else v) for k, v in row.items()}
                for i, row in allper.items()}

    vi = slice_arm(0, ntre)
    ri = slice_arm(ntre, 2 * ntre)
    ci = slice_arm(2 * ntre, 2 * ntre + nctrl)

    published = {}
    for r in csv.DictReader((analysis / 'lagrangian_budget.csv').open()):
        i = int(r['index'])
        if SEED_INDEX < i <= LAST_INDEX:
            published.setdefault(i, {})[int(r['parcel'])] = r
    diffs = {}
    for i, row in vi.items():
        for p in range(ntre):
            for k in keys:
                if k not in published[i][p]:
                    continue
                a = row[k][p]
                b = float(published[i][p][k])
                if np.isfinite(a) and np.isfinite(b):
                    diffs[k] = max(diffs.get(k, 0.0), abs(a - b))
    summary['sampler_validation_max_abs_diff'] = {k: float(v) for k, v in diffs.items()}
    summary['sampler_validation_worst'] = float(max(diffs.values())) if diffs else None
    print('  max abs diff vs published: %r' % summary['sampler_validation_worst'], flush=True)

    track = {int(r['index']): r for r in csv.DictReader((analysis / 'vortex_track.csv').open())
             if r.get('zeta_low_max_s')}
    peak = track[93]
    vx, vy = float(peak['center_x_m']), float(peak['center_y_m'])
    zeta93 = curl(native(groups[93], NK), coords)[2]

    def zsample(a, p):
        return RegularGridInterpolator(coords, a, bounds_error=False, fill_value=np.nan)(p)

    rows = []
    for arm, pos, per in (('treatment', treatment, vi), ('retrace', retrace_pos, ri),
                          ('control', ctrl_pos, ci)):
        n = len(pos[SEED_INDEX])
        pi = phase_integrals(per, keys, n)
        p93 = pos[93]
        d = np.hypot(p93[:, 0] - vx, p93[:, 1] - vy)
        z93 = zsample(zeta93, p93)
        for p in range(n):
            seed = pos[SEED_INDEX][p]
            base = dict(arm=arm, parcel=p,
                        seed_x_m=float(seed[0]), seed_y_m=float(seed[1]), seed_z_m=float(seed[2]),
                        stratum=strata[p][0] if arm == 'control' else arm,
                        seed_radius_m=(strata[p][1] if arm == 'control'
                                       else float(np.hypot(seed[0] - centre[0],
                                                           seed[1] - centre[1]))),
                        dist_to_vortex_2790_m=float(d[p]), zeta_2790_s=float(z93[p]),
                        z_2790_m=float(p93[p, 2]))
            for a, b, acc, valid in pi:
                row = dict(base)
                row.update(phase='%.0f-%.0f' % (times[a], times[b]),
                           phase_start_s=float(times[a]), phase_end_s=float(times[b]),
                           valid=bool(valid[p]))
                for k in keys:
                    row[k] = float(acc[k][p])
                rows.append(row)
    write_csv(out / 'control_phase_integrals.csv', rows)
    summary['rows'] = len(rows)
    summary['control_n'] = len(seeds)
    summary['operator_keys'] = keys
    (out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2), flush=True)
    f.close()


if __name__ == '__main__':
    main()
