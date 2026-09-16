"""Summarize the completed discrete audit without rereading simulation fields."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT/'outputs/les_discrete_balance_20260915_v2')
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/les_discrete_balance_synthesis_20260915')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    source_files = [args.input / name for name in
                    ('terms_by_height.csv', 'terms_by_block.csv', 'balances_by_block.csv', 'summary.json')]
    totals, blocks, balances = [read_csv(p) for p in source_files[:3]]
    source = json.loads(source_files[3].read_text())
    level = 121.65909956343168
    select = lambda rows, dx, radius, mode, quantity: [r for r in rows
        if int(r['grid_dx_m']) == dx and float(r['radius_m']) == radius
        and r['mask_mode'] == mode and r['quantity'] == quantity
        and abs(float(r['level_m']) - level) < 1e-8]
    cases = {}
    radial = []
    for dx in (600, 300):
        modes = {}
        for mode in ('moving_end', 'moving_symmetric', 'fixed_initial'):
            modes[mode] = {}
            for quantity in ('les_label', 'total'):
                rows = select(totals, dx, 4200, mode, quantity)
                values = {r['term']:float(r['signed_m2_s']) for r in rows}
                values['sum_m2_s'] = sum(values.values())
                modes[mode][quantity] = values
        for radius in (1200,2400,4200,6000):
            rows = sorted(select(balances,dx,radius,'moving_end','total'),key=lambda r:int(r['block']))
            radial.append(dict(grid_dx_m=dx,radius_m=radius,
                initial_m2_s=float(rows[0]['inventory_start_m2_s']),
                final_m2_s=float(rows[-1]['inventory_end_m2_s']),
                change_m2_s=float(rows[-1]['inventory_end_m2_s'])-float(rows[0]['inventory_start_m2_s'])))
        by_radius = {r['radius_m']:r for r in radial if r['grid_dx_m']==dx}
        annulus = {k: by_radius[4200][k]-by_radius[1200][k]
                   for k in ('initial_m2_s','final_m2_s','change_m2_s')}
        op_signs = {}
        for mode in modes:
            subset = select(blocks,dx,4200,mode,'total')
            op_signs[mode] = {term:sum(float(r['signed_m2_s'])<0 for r in subset if r['term']==term)
                              for term in ('advection','les','mask_motion')}
        context = source['cases'][str(dx)]['context_from_existing_tracking']
        zpeak = context['zeta_max_s-1']['maximum_time_s']
        cpeak = context['convergence_max_s-1']['maximum_time_s']
        first = source['cases'][str(dx)]['blocks'][0]['time_start_s']
        cases[str(dx)] = dict(level_m=level, modes_at_4200m=modes,
            annulus_1200_to_4200_m=annulus, negative_blocks_of_nine=op_signs,
            peak_convergence_minus_peak_zeta_s=cpeak-zpeak,
            peak_at_left_boundary=(zpeak==first or cpeak==first),
            context=context)

    fig, axes = plt.subplots(1,3,figsize=(15,4.8),constrained_layout=True)
    colors = {600:'#167b72',300:'#b43850'}
    labels = ['MUSCL total','LES local','Coriolis','Arrasto','Mascara']
    keys = ['advection','les','coriolis','surface_drag','mask_motion']
    for dx, offset in ((600,-.18),(300,.18)):
        q = cases[str(dx)]['modes_at_4200m']['moving_end']
        axes[0].bar(np.arange(5)+offset,[q['total'][k]/1000 for k in keys],
                    width=.34,color=colors[dx],label=f'{dx} m')
        lk = ['local_les','inferred_muscl_label','mask_motion']
        axes[1].bar(np.arange(3)+offset,[q['les_label'][k]/1000 for k in lk],width=.34,color=colors[dx])
        data = [r for r in radial if r['grid_dx_m']==dx]
        axes[2].plot([r['radius_m']/1000 for r in data],[r['initial_m2_s']/1000 for r in data],
                     '--o',color=colors[dx],label=f'{dx} m inicial')
        axes[2].plot([r['radius_m']/1000 for r in data],[r['final_m2_s']/1000 for r in data],
                     '-s',color=colors[dx],label=f'{dx} m final')
    axes[0].set_xticks(range(5),labels,rotation=25,ha='right')
    axes[1].set_xticks(range(3),['LES local','Evolucao inferida\ndo rotulo','Mascara'])
    axes[0].set_title('Circulacao total: termos acumulados')
    axes[1].set_title('Inventario rotulado LES: termos acumulados')
    axes[2].set(title='Circulacao assinada por raio',xlabel='Raio (km)')
    for ax in axes:
        ax.axhline(0,color='#333333',linewidth=.7)
        ax.set_ylabel('10^3 m2/s')
        ax.grid(axis='y',alpha=.2)
        ax.set_axisbelow(True)
    axes[0].legend()
    axes[2].legend(fontsize=8)
    fig.suptitle('Fase madura, z = 121,7 m | disco movel, convencao de ponto final | r = 4,2 km nos balancos')
    fig.savefig(args.out/'circulation_balance.png',dpi=170)
    plt.close(fig)
    with (args.out/'radial_circulation.csv').open('x',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(radial[0]));writer.writeheader();writer.writerows(radial)
    result=dict(cases=cases, source_summary=str(source_files[3].relative_to(ROOT)),
                interpretation='operator accounting on saved trajectories, not physical intervention')
    (args.out/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(input_sha256={str(p.relative_to(ROOT)):digest(p) for p in source_files},
                  script_sha256=digest(Path(__file__)),
                  artifacts={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in args.out.iterdir()})
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({dx:dict(annulus=c['annulus_1200_to_4200_m'],
          peak_offset_s=c['peak_convergence_minus_peak_zeta_s']) for dx,c in cases.items()},indent=2))


if __name__ == '__main__':
    main()
