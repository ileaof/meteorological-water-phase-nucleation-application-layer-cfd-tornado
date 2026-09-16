"""Audit temporal cancellation in saved LES inventory balances; no simulation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from statistics import median


def audit(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(float(row['grid_dx_m']), int(row['level_index']))].append(row)
    result = []
    for (dx, level), group in sorted(groups.items()):
        production = [float(r['local_production_signed']) for r in group]
        residual = [float(r['transport']) for r in group]
        gross = sum(float(r['local_production_abs']) for r in group)
        p_net = sum(production)
        r_net = sum(residual)
        p_abs = sum(map(abs, production))
        r_abs = sum(map(abs, residual))
        error = max(abs(float(r['inventory_change']) - p - t)
                    for r, p, t in zip(group, production, residual))
        result.append(dict(
            grid_dx_m=dx, level_index=level, level_m=float(group[0]['level_m']),
            intervals=len(group), production_net=p_net, production_abs_intervals=p_abs,
            production_spatial_abs_intervals=gross, residual_net=r_net,
            residual_abs_intervals=r_abs,
            old_net_residual_share=abs(r_net) / (gross + abs(r_net)),
            uncancelled_residual_share=r_abs / (gross + r_abs),
            signed_balance_residual_share=r_abs / (p_abs + r_abs),
            residual_retained_fraction=abs(r_net) / r_abs if r_abs else 0.0,
            balance_identity_max_abs=error,
        ))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(
        'outputs/les_local_vs_transported_20260913/les_local_vs_transported.csv'))
    parser.add_argument('--out', type=Path, default=Path(
        'outputs/les_residual_cancellation_20260915'))
    args = parser.parse_args()
    with args.input.open(newline='', encoding='utf-8') as stream:
        rows = audit(list(csv.DictReader(stream)))
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / 'by_height.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metrics = ('old_net_residual_share', 'uncancelled_residual_share',
               'signed_balance_residual_share', 'residual_retained_fraction')
    summary = {'input': str(args.input),
               'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
               'by_resolution': {},
               'limitations': [
                   'Residual is inventory change minus local LES increment; not a measured boundary flux.',
                   'Absolute values are taken after spatial integration; spatial cancellation remains.',
                   'Gross production uses interval-integrated increments; native-step cancellation remains.',
                   'Moving masks are fixed within each difference, but change between intervals.',
                   '600 m source intervals are not exactly aligned; no error bound follows from time offset alone.',
                   'Balance identity is arithmetic consistency, not independent physical validation.',
               ]}
    for dx in sorted({r['grid_dx_m'] for r in rows}):
        subset = [r for r in rows if r['grid_dx_m'] == dx]
        summary['by_resolution'][str(int(dx))] = {
            metric: {'min': min(r[metric] for r in subset),
                     'median': median(r[metric] for r in subset),
                     'max': max(r[metric] for r in subset)} for metric in metrics}
    summary['balance_identity_max_abs'] = max(r['balance_identity_max_abs'] for r in rows)
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
