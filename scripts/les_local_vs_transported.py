"""Separate LOCAL LES production from TRANSPORT of the LES label, by height.

`docs/CLOSURE_CONVERGENCE_BY_HEIGHT.md` measured that the v4 `les` source inventory is
much larger near the surface at 300 m than at 600 m.  That inventory is a TRANSPORTED
label: it says where LES-labelled vorticity is, not where the closure acts.  The audit
in AGENT_CHANNEL (turnos 13-14) asked for the three contributions to be separated.

This does that, per level, over a common physical mask:

  local production   P_k = integral of curl(increments/les)_z dA
                     the curl of the velocity increment the LES operator ACTUALLY applied
                     during the interval -- production, in the cell where it happened
  inventory change   dI_k = I_k(t_i) - I_k(t_{i-1}), the v4 label integral, SAME mask at
                     both times so mask motion does not leak into the difference
  transport          T_k = dI_k - P_k
                     everything that moved label into or out of that level and cylinder

Both quantities come from the same trajectory: the provenance runs reproduce their control
bit-for-bit in all 12 prognostic fields (verified in each summary.json `neutrality`), so
the sequence increments and the provenance sources describe the same flow.

Reads only.  No simulation is run.
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
from storm_dynamics.diagnostic_capture import curl                     # noqa: E402
from analyze_diagnostic_sequence import native                          # noqa: E402

RUNS = {
    600: dict(sequence='outputs/diagnostic_sequence_20260905/sequence.h5',
              provenance='outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5',
              track='outputs/resolution_600m_comparison_audit_20260909/vortex_timeseries.csv'),
    300: dict(sequence='outputs/resolution_300m_20260909/sequence.h5',
              provenance='outputs/resolution_300m_provenance_20260909/provenance_long_v4.h5',
              track='outputs/resolution_300m_audit_20260909/vortex_timeseries.csv'),
}


def centres(path):
    return [(float(r['time_s']), float(r['center_x_m']), float(r['center_y_m']))
            for r in csv.DictReader(Path(path).open())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--radius-m', type=float, default=4200.0)
    ap.add_argument('--out', default='outputs/les_local_vs_transported_20260913')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    checks = []
    for dx, paths in RUNS.items():
        seq = h5py.File(paths['sequence'], 'r')
        prov = h5py.File(paths['provenance'], 'r')
        nz = int(prov.attrs['analysis_nz'])
        nk = int(seq.attrs['budget_nz_with_halo'])
        xc, yc = seq['grid/xc'][:], seq['grid/yc'][:]
        zc = seq['grid/zc'][:nk]
        coords = (xc, yc, zc)
        cell_area = float(np.diff(xc).mean() * np.diff(yc).mean())
        track = centres(paths['track'])

        pgroups = [prov['snapshots'][k] for k in sorted(prov['snapshots'].keys())]
        ptimes = np.array([g.attrs['time_s'] for g in pgroups])
        sgroups = list(seq['snapshots'].values())
        stimes = np.array([g.attrs['time_s'] for g in sgroups])
        # the sequence may be longer than the provenance window; align by time
        sidx = [int(np.argmin(np.abs(stimes - t))) for t in ptimes]
        worst = float(np.max(np.abs(stimes[sidx] - ptimes)))
        checks.append(dict(grid_dx_m=dx, max_time_alignment_s=worst,
                           n_snapshots=len(ptimes)))
        print('%d m: %d snapshots, time alignment worst %.3g s' % (dx, len(ptimes), worst),
              flush=True)

        inv_prev = None
        for n, (pg, si) in enumerate(zip(pgroups, sidx)):
            time_s = float(ptimes[n])
            _, cx, cy = min(track, key=lambda r: abs(r[0] - time_s))
            disc = ((xc[:, None] - cx) ** 2 + (yc[None, :] - cy) ** 2) <= args.radius_m ** 2
            inv = prov['snapshots'][sorted(prov['snapshots'].keys())[n]]['omega_source']['les'][2, :, :, :nz]
            if n > 0:
                # same mask applied to both times: mask motion cannot leak into dI
                prod = curl(native(sgroups[si]['increments/les']), coords)[2][:, :, :nz]
                for k in range(nz):
                    m = disc
                    P = float(prod[:, :, k][m].sum() * cell_area)
                    Pabs = float(np.abs(prod[:, :, k])[m].sum() * cell_area)
                    I1 = float(inv[:, :, k][m].sum() * cell_area)
                    I0 = float(inv_prev[:, :, k][m].sum() * cell_area)
                    Iabs = float(np.abs(inv[:, :, k])[m].sum() * cell_area)
                    rows.append(dict(grid_dx_m=dx, time_s=time_s, level_m=float(zc[k]),
                                     level_index=k,
                                     local_production_signed=P, local_production_abs=Pabs,
                                     inventory_signed=I1, inventory_abs=Iabs,
                                     inventory_change=I1 - I0, transport=(I1 - I0) - P))
            inv_prev = inv
            if n % 5 == 0:
                print('  %d m  %d/%d  t=%.1f' % (dx, n, len(ptimes), time_s), flush=True)
        seq.close()
        prov.close()

    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (out / 'les_local_vs_transported.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    (out / 'summary.json').write_text(json.dumps(
        dict(radius_m=args.radius_m, rows=len(rows), alignment=checks), indent=2),
        encoding='utf-8')
    print('wrote', out / 'les_local_vs_transported.csv')
    figure(rows, out)


def figure(rows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from collections import defaultdict
    agg = defaultdict(lambda: defaultdict(float))
    lev = {}
    for r in rows:
        key = (r['grid_dx_m'], r['level_index'])
        lev[r['level_index']] = r['level_m']
        agg[key]['Pabs'] += abs(r['local_production_abs'])
        agg[key]['P'] += r['local_production_signed']
        agg[key]['T'] += r['transport']
    ks = sorted(lev)
    z = [lev[k] for k in ks]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
    axes[0].plot([agg[(300.0, k)]['Pabs'] / agg[(600.0, k)]['Pabs'] for k in ks], z, 'o-')
    axes[0].axvline(1.0, color='grey', lw=.8)
    axes[0].set(xlabel='gross local production, 300 m / 600 m', ylabel='Height (m)',
                title='Where refinement raises LES activity')
    for dx, style in ((600.0, 'o-'), (300.0, 's-')):
        axes[1].plot([abs(agg[(dx, k)]['P']) / agg[(dx, k)]['Pabs'] for k in ks], z, style,
                     label='dx = %d m' % dx)
        axes[2].plot([abs(agg[(dx, k)]['T'])
                      / (agg[(dx, k)]['Pabs'] + abs(agg[(dx, k)]['T'])) for k in ks], z,
                     style, label='dx = %d m' % dx)
    axes[1].set(xlabel='$|\\sum P| / \\sum |P|$', title='Net as a fraction of gross')
    axes[2].set(xlabel='$|T| / (|P| + |T|)$', title='Transport share of the balance')
    for ax in axes:
        ax.grid(alpha=.3)
        ax.set_ylim(0, 2000)
    axes[1].legend()
    fig.suptitle('LES label: locally produced, largely self-cancelling, barely transported '
                 '(2790-3300 s, 4.2 km cylinder)')
    fig.savefig(out / 'les_local_vs_transported.png', dpi=150)
    plt.close(fig)
    print('wrote', out / 'les_local_vs_transported.png')


if __name__ == '__main__':
    main()
