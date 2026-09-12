"""Stratified table and figure for the control-parcel case-control test.

Reads control_phase_integrals.csv produced by control_parcel_test.py and writes
the outcome-stratified comparison used in docs/CONTROL_PARCEL_TEST.md.

The height bands are a DESCRIPTIVE stratifier applied after the fact; the
seeding rule in control_parcel_test.py stayed blind to the outcome.
"""
from pathlib import Path
import argparse
import csv

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

BANDS = [(0, 400), (400, 700), (700, 1200), (1200, 3000)]
TERMS = ('les_zeta', 'tilting_zeta', 'stretching_zeta')


def load(path):
    rows = list(csv.DictReader(Path(path).open()))
    for r in rows:
        for k, v in r.items():
            if k in ('arm', 'stratum', 'phase'):
                continue
            r[k] = (v == 'True') if k == 'valid' else float(v)
    return rows


def describe(sel):
    out = {'n': len(sel)}
    if not sel:
        return out
    arr = {k: np.array([r[k] for r in sel]) for k in TERMS + ('delta_zeta',)}
    for k, v in arr.items():
        out[k + '_median'] = float(np.median(v))
        out[k + '_positive'] = int((v > 0).sum())
    mag = sum(np.abs(arr[k]) for k in TERMS)
    out['les_share'] = float(np.median(np.abs(arr['les_zeta']) / mag))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='outputs/control_parcels_20260912_dense/'
                                       'control_phase_integrals.csv')
    ap.add_argument('--out', default='outputs/control_parcels_20260912_dense')
    args = ap.parse_args()
    out = Path(args.out)
    rows = load(args.input)
    phases = sorted({r['phase'] for r in rows}, key=lambda p: float(p.split('-')[0]))

    table = []
    for phase in phases:
        pick = [r for r in rows if r['phase'] == phase and r['valid']]
        groups = [('treatment', [r for r in pick if r['arm'] == 'treatment'])]
        vortex = [r for r in pick if r['arm'] == 'control' and r['zeta_2790_s'] >= 0.003]
        for lo, hi in BANDS:
            groups.append(('control %d-%d m' % (lo, hi),
                           [r for r in vortex if lo <= r['z_2790_m'] < hi]))
        groups.append(('control non-vortex',
                       [r for r in pick if r['arm'] == 'control' and r['zeta_2790_s'] < 0.003]))
        for name, sel in groups:
            table.append(dict(phase=phase, group=name, **describe(sel)))
    keys = list(dict.fromkeys(k for row in table for k in row))
    with (out / 'stratified_comparison.csv').open('w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(table)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), constrained_layout=True)
    centres = [(lo + hi) / 2 for lo, hi in BANDS]
    for ax, phase, title in zip(axes, (phases[1], phases[2]),
                                ('Preparatory phase 2370-2610 s',
                                 'Amplification phase 2610-2790 s')):
        sub = [r for r in table if r['phase'] == phase and r['group'].startswith('control ')
               and r['group'] != 'control non-vortex']
        for key, label in (('les_zeta', 'LES closure'), ('tilting_zeta', 'tilting'),
                           ('stretching_zeta', 'stretching')):
            ax.plot(centres, [r[key + '_median'] for r in sub], 'o-', label=label)
        tre = [r for r in table if r['phase'] == phase and r['group'] == 'treatment'][0]
        for key, colour in (('les_zeta', 'C0'), ('tilting_zeta', 'C1'), ('stretching_zeta', 'C2')):
            ax.plot([492], [tre[key + '_median']], marker='*', ms=16, color=colour,
                    linestyle='none')
        ax.axhline(0, color='grey', lw=.8)
        ax.set(xlabel='Parcel height at 2790 s (m)',
               ylabel='Phase integral of $\\zeta$ (s$^{-1}$)', title=title)
        ax.grid(alpha=.3)
    axes[0].legend()
    fig.suptitle('Blind-seeded control parcels that reach the vortex; stars = the 27 '
                 'outcome-selected parcels')
    fig.savefig(out / 'control_stratified.png', dpi=150)
    plt.close(fig)
    print('wrote', out / 'stratified_comparison.csv', 'and', out / 'control_stratified.png')


if __name__ == '__main__':
    main()
