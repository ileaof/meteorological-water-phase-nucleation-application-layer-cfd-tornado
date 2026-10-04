"""Analyze a gated direct capture; never advances a simulation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def write_csv(path,rows):
    with path.open('x',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=ROOT/'outputs/muscl_direct_capture_20260920')
    parser.add_argument('--out',type=Path,default=ROOT/'outputs/muscl_direct_analysis_20260920')
    args=parser.parse_args()
    meta=json.loads((args.input/'metadata.json').read_text())
    if meta['status']!='complete' or meta['completed_steps']!=1015 or not meta['pilot_pass']:
        raise RuntimeError('Complete gated replay required for mature-window conclusions')
    args.out.mkdir(parents=True,exist_ok=False)
    with h5py.File(args.input/'capture.h5','r') as f:
        assert f.attrs['status']=='complete' and int(f.attrs['completed_blocks'])==9
        terms=json.loads(f.attrs['terms']);masks=json.loads(f.attrs['masks'])
        steps=f['step_statistics'][:];blocks=f['block_statistics'][:]
        balances=f['block_balances'][:];levels=f['levels_m'][:];schedule=f['schedule'][:]
        block_steps=[(int(g.attrs['start_step']),int(g.attrs['end_step'])) for g in f['blocks'].values()]
    assert np.isfinite(steps).all() and np.isfinite(blocks).all() and np.isfinite(balances).all()
    rows=[];blockrows=[]
    for t,term in enumerate(terms):
        for m,mask in enumerate(masks):
            for k,level in enumerate(levels):
                signed=steps[:,t,m,k,0]
                spatial=steps[:,t,m,k,1]
                block_signed=blocks[:,t,m,k,0]
                block_spatial=blocks[:,t,m,k,1]
                np.testing.assert_allclose(signed.sum(),block_signed.sum(),rtol=1e-10,atol=1e-8)
                rows.append(dict(term=term,**mask,level_m=float(level),
                    signed_sum_m2_s=float(signed.sum()),
                    sum_abs_step_integrals_m2_s=float(np.abs(signed).sum()),
                    sum_step_spatial_abs_m2_s=float(spatial.sum()),
                    sum_abs_block_integrals_m2_s=float(np.abs(block_signed).sum()),
                    sum_block_spatial_abs_m2_s=float(block_spatial.sum()),
                    negative_steps=int((signed<0).sum()),negative_blocks=int((block_signed<0).sum()),
                    spatial_retained_at_step_cadence=float(np.abs(signed).sum()/spatial.sum()) if spatial.sum() else None,
                    temporal_retained_signed=float(abs(signed.sum())/np.abs(signed).sum()) if np.abs(signed).sum() else None,
                    block_vs_step_spatial_retained=float(block_spatial.sum()/spatial.sum()) if spatial.sum() else None))
                for b in range(9):
                    blockrows.append(dict(block=b,term=term,**mask,level_m=float(level),
                                          signed_m2_s=float(block_signed[b]),spatial_abs_m2_s=float(block_spatial[b])))
    write_csv(args.out/'directional_and_cancellation.csv',rows)
    write_csv(args.out/'terms_by_block.csv',blockrows)
    k=int(np.argmin(abs(levels-121.65909956343168)))
    primary=[r for r in rows if r['level_m']==float(levels[k])]
    write_csv(args.out/'low_level_all_radii.csv',primary)
    annuli=[]
    for mode in ('moving_end','moving_symmetric','fixed_initial'):
        for inner,outer in ((1200,4200),(4200,6000)):
            mi=next(i for i,m in enumerate(masks) if m['radius_m']==inner and m['mask_mode']==mode)
            mo=next(i for i,m in enumerate(masks) if m['radius_m']==outer and m['mask_mode']==mode)
            for t,term in enumerate(terms):
                for j,level in enumerate(levels):
                    signed=steps[:,t,mo,j,0]-steps[:,t,mi,j,0]
                    spatial=steps[:,t,mo,j,1]-steps[:,t,mi,j,1]
                    annuli.append(dict(term=term,mask_mode=mode,inner_radius_m=inner,
                        outer_radius_m=outer,level_m=float(level),signed_sum_m2_s=float(signed.sum()),
                        sum_abs_step_integrals_m2_s=float(abs(signed).sum()),
                        sum_step_spatial_abs_m2_s=float(spatial.sum()),
                        negative_steps=int((signed<0).sum())))
    write_csv(args.out/'annular_contributions.csv',annuli)
    inventory=[]
    for m,mask in enumerate(masks):
        for q,quantity in enumerate(('total','les_label')):
            for j,level in enumerate(levels):
                inventory.append(dict(**mask,quantity=quantity,level_m=float(level),
                    inventory_start_m2_s=float(balances[0,m,q,j,0]),
                    inventory_end_m2_s=float(balances[-1,m,q,j,1]),
                    delta_inventory_m2_s=float((balances[:,m,q,j,1]-balances[:,m,q,j,0]).sum()),
                    mask_motion_m2_s=float(balances[:,m,q,j,2].sum())))
    write_csv(args.out/'inventory_and_mask_motion.csv',inventory)
    counts={n:sum(all(v['bitwise'] for v in s['fields'].values()) for s in meta['original_snapshots'])
            for n in ('original_snapshots_all_12_bitwise',)}
    summary=dict(status='complete',level_m=float(levels[k]),steps=1015,blocks=9,
        gates=meta['gates'],**counts,neutrality_bitwise_all_12=meta['neutrality_bitwise_all_12'],
        old_label_checks=meta['old_label_checks'],wall_clock_s=meta['wall_clock_s'],
        maximum_prior_balance_difference_m2_s=max(r['max_abs_m2_s'] for r in meta['prior_balance_comparison']),
        low_level_results=primary,
        limitations=[
            'Directional terms are momentum-flux contributions, not physical vorticity fluxes.',
            'Flux formulas are reevaluated by the native routines; actual updates provide a separate reconstruction check.',
            'Masks are fixed within each archived block; no native-step moving-centre trajectory was invented.',
            'Other total operators are combined within each step, retaining inter-stage cancellation.',
            'Nine blocks are not independent realizations; magnitudes are not causal percentages.',
            'Only the mature 600 m trajectory is replayed; no convergence or counterfactual physics experiment.'])
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    times=schedule[:,0]+schedule[:,1]
    fig,axes=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    curves=(('total','moving_end'),('total','moving_symmetric'),('total','fixed_initial'),('les','moving_end'))
    for ax,(prefix,mode) in zip(axes.ravel(),curves):
        m=next(i for i,x in enumerate(masks) if x['radius_m']==4200 and x['mask_mode']==mode)
        for direction,color in zip('xyz',('#167b72','#ba465c','#926a16')):
            t=terms.index(prefix+'_flux_'+direction)
            ax.plot(times,np.cumsum(steps[:,t,m,k,0])/1000,label=direction,color=color)
        t=terms.index(prefix+'_muscl')
        ax.plot(times,np.cumsum(steps[:,t,m,k,0])/1000,label='soma medida',color='#222222',ls='--')
        ax.set(title=f'{prefix} | {mode}',xlabel='Tempo simulado (s)',ylabel='Incremento acumulado (10^3 m2/s)')
        ax.axhline(0,color='gray',lw=.6);ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Fluxos de momento por direcao | 600 m | z=121,7 m | raio 4,2 km')
    fig.savefig(args.out/'directional_circulation.png',dpi=160);plt.close(fig)
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    manifest=dict(input_sha256={str(p.relative_to(ROOT)):digest(p) for p in
        (args.input/'metadata.json',args.input/'capture.h5')},script_sha256=digest(Path(__file__)),
        artifacts={p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in args.out.iterdir()})
    (args.out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps([r for r in primary if r['radius_m']==4200 and r['mask_mode']=='moving_end'],indent=2))


if __name__=='__main__':main()
