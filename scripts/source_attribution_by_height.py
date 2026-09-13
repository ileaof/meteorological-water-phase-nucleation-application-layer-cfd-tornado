"""Height-stratified vorticity source attribution, 600 m vs 300 m.

The published provenance integrals aggregate each source over a single vortex
mask.  The control-parcel test (docs/CONTROL_PARCEL_TEST.md) showed that height
is the stratifier that decides whether the LES closure or the resolved terms
dominate, so the aggregate hides the quantity of interest.

This recomputes the v4 additive source attribution LEVEL BY LEVEL, on the
matched 600 m / 300 m pair that covers the same window from the same state, so
the near-surface LES share can be compared across a factor of two in dx.

Both runs share an identical vertical grid; only dx changes.  The mask is a
cylinder of fixed PHYSICAL radius centred on each run's own tracked vortex, so
the comparison is in metres, never in cells.

Reads only.  No simulation is run.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
import h5py

RUNS = {
    600: ('outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5',
          'outputs/resolution_600m_comparison_audit_20260909/vortex_timeseries.csv'),
    300: ('outputs/resolution_300m_provenance_20260909/provenance_long_v4.h5',
          'outputs/resolution_300m_audit_20260909/vortex_timeseries.csv'),
}


def centres(path):
    out = []
    for r in csv.DictReader(Path(path).open()):
        out.append((float(r['time_s']), float(r['center_x_m']), float(r['center_y_m'])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--radius-m', type=float, default=4200.0)
    ap.add_argument('--out', default='outputs/source_attribution_by_height_20260913')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    closure = []
    for dx, (h5path, csvpath) in RUNS.items():
        f = h5py.File(h5path, 'r')
        names = json.loads(f.attrs['source_names'])
        nz = int(f.attrs['analysis_nz'])
        xc = f['grid/xc'][:]
        yc = f['grid/yc'][:]
        zc = f['grid/zc'][:]
        cell_area = float(np.diff(xc).mean() * np.diff(yc).mean())
        track = centres(csvpath)
        keys = sorted(f['snapshots'].keys())
        print('%d m: %d snapshots, %d sources, %d levels' % (dx, len(keys), len(names), nz),
              flush=True)
        for n, key in enumerate(keys):
            g = f['snapshots'][key]
            time_s = float(g.attrs['time_s'])
            _, cx, cy = min(track, key=lambda r: abs(r[0] - time_s))
            r2 = ((xc[:, None] - cx) ** 2 + (yc[None, :] - cy) ** 2)
            disc = r2 <= args.radius_m ** 2
            total = g['omega_total'][2, :, :, :nz]
            stack = {}
            for s in names:
                stack[s] = g['omega_source'][s][2, :, :, :nz]
            summed = sum(stack.values())
            scale = max(float(np.abs(total).max()), 1e-30)
            closure.append(dict(grid_dx_m=dx, time_s=time_s,
                                partition_max_abs_error=float(np.abs(summed - total).max()),
                                partition_relative=float(np.abs(summed - total).max() / scale)))
            for k in range(nz):
                m = disc
                tot_abs = float(np.abs(total[:, :, k])[m].sum() * cell_area)
                tot_sig = float(total[:, :, k][m].sum() * cell_area)
                denom = 0.0
                per = {}
                for s in names:
                    a = float(np.abs(stack[s][:, :, k])[m].sum() * cell_area)
                    per[s] = (float(stack[s][:, :, k][m].sum() * cell_area), a)
                    denom += a
                for s in names:
                    sig, ab = per[s]
                    rows.append(dict(grid_dx_m=dx, time_s=time_s, level_m=float(zc[k]),
                                     level_index=k, source=s,
                                     signed_zeta_integral=sig, absolute_zeta_integral=ab,
                                     absolute_share=ab / denom if denom else np.nan,
                                     total_signed=tot_sig, total_absolute=tot_abs))
            if n % 5 == 0:
                print('  %d m  %d/%d  t=%.1f' % (dx, n, len(keys), time_s), flush=True)
        f.close()

    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (out / 'source_by_height.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    with (out / 'partition_closure.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(closure[0]))
        w.writeheader()
        w.writerows(closure)
    worst = max(c['partition_relative'] for c in closure)
    (out / 'summary.json').write_text(json.dumps(
        dict(radius_m=args.radius_m, rows=len(rows),
             partition_relative_max=worst,
             runs={str(k): v[0] for k, v in RUNS.items()}), indent=2), encoding='utf-8')
    print('partition closure, worst relative error:', worst)
    print('wrote', out / 'source_by_height.csv')
    figure(rows, out)


def figure(rows, out):
    """One profile per resolution at ITS OWN final snapshot.

    An earlier version selected every output within 40 s of the final time.  The
    snapshot spacing is 30 s, so that took TWO instants per height and joined
    them as if they were one profile.  Select the exact final time per run.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    series = {}
    finals = {}
    for dx in (600, 300):
        finals[dx] = max(r['time_s'] for r in rows if r['grid_dx_m'] == dx)
        pick = [r for r in rows if r['grid_dx_m'] == dx and r['source'] == 'les'
                and r['time_s'] == finals[dx]]
        pick.sort(key=lambda r: r['level_index'])
        series[dx] = ([r['level_m'] for r in pick], [r['absolute_share'] for r in pick],
                      [r['absolute_zeta_integral'] for r in pick])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    for dx, style in ((600, 'o-'), (300, 's-')):
        z, sh, ab = series[dx]
        axes[0].plot(sh, z, style, label='dx = %d m' % dx)
        axes[1].plot(ab, z, style, label='dx = %d m' % dx)
    ratio = [b / a if a else np.nan for a, b in zip(series[600][1], series[300][1])]
    axes[2].plot(ratio, series[600][0], 'd-', color='C3')
    axes[2].axvline(1.0, color='grey', lw=.8)
    axes[0].set(xlabel='LES share of $|\\zeta|$ sources', ylabel='Height (m)',
                title='LES-labelled share by height')
    axes[1].set(xlabel='$\\int|\\zeta_{LES}|\\,dA$ (m$^2$ s$^{-1}$)',
                title='Absolute LES-labelled $\\zeta$')
    axes[2].set(xlabel='share at 300 m / share at 600 m', title='Effect of halving dx')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.set_ylim(0, 2000)
    axes[0].legend()
    fig.suptitle('LES-labelled vorticity INVENTORY (transported since restart), not local '
                 'production — t = %.1f s, 4.2 km cylinder' % finals[300])
    fig.savefig(out / 'closure_share_by_height.png', dpi=150)
    plt.close(fig)
    print('wrote', out / 'closure_share_by_height.png')


if __name__ == '__main__':
    main()
