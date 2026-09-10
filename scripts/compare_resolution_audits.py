"""Compare gated audits on common physical sample times without extrapolation."""
from pathlib import Path
import csv
import json
import ast
import hashlib
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/resolution_comparison_20260909'
CASES={600:ROOT/'outputs/resolution_600m_comparison_audit_20260909',
       300:ROOT/'outputs/resolution_300m_audit_20260909'}


def read(path):
    with path.open(encoding='utf-8') as f: return list(csv.DictReader(f))


def write(path,rows):
    if not rows: return
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    summaries={dx:json.loads((p/'summary.json').read_text()) for dx,p in CASES.items()}
    for s in summaries.values():
        assert s['numerical_gate']=='PASS'
        assert s['native_budget']['max_closure_relative_rms']<1e-10
    expected_geometry={'radial_bin_m':600.0,'cylinder_radius_m':4200.0,'vertical_range_m':[0.0,2000.0]}
    for dx,s in summaries.items():
        assert s.get('analysis_geometry')==expected_geometry,(dx,s.get('analysis_geometry'))
    data={dx:read(p/'vortex_timeseries.csv') for dx,p in CASES.items()}
    ts={dx:np.array([float(r['time_s']) for r in rows]) for dx,rows in data.items()}
    first=max(t[0] for t in ts.values());last=min(t[-1] for t in ts.values())
    common=np.r_[np.arange(np.ceil(first/30)*30,last-1e-9,30),last]
    common=np.unique(np.round(common,9))
    excluded={'time_s','sequence_time_s','alignment_time_mismatch_s','center_x_m','center_y_m'}
    shared=set(data[600][0]) & set(data[300][0])
    metrics=sorted(k for k in shared-excluded if not k.endswith('_cells'))
    cell_metrics=sorted(k for k in shared-excluded if k.endswith('_cells'))
    paired=[];stats=[];adequacy=[]
    for key in metrics+cell_metrics:
        values={dx:np.interp(common,ts[dx],[float(r[key]) for r in rows]) for dx,rows in data.items()}
        a,b=values[600],values[300]
        valid=np.isfinite(a)&np.isfinite(b)
        target=adequacy if key.endswith('_cells') else stats
        target.append(dict(metric=key,median_600m=float(np.nanmedian(a)),median_300m=float(np.nanmedian(b)),
                           median_paired_delta=float(np.nanmedian(b-a)),final_600m=float(a[-1]),final_300m=float(b[-1]),
                           fraction_300m_greater=float(np.mean(b[valid]>a[valid])) if valid.any() else float('nan'),
                           valid_pairs=int(valid.sum())))
        if not key.endswith('_cells'):
            for t,x,y in zip(common,a,b):
                paired.append(dict(time_s=float(t),metric=key,value_600m=float(x),value_300m=float(y),delta=float(y-x)))
    write(OUT/'paired_timeseries.csv',paired);write(OUT/'metric_comparison.csv',stats)
    write(OUT/'resolution_adequacy_cells.csv',adequacy)
    # Preserve the full source and budget definitions alongside the scalar comparison.
    for filename,groupkeys in [('provenance_integrals.csv',('source',)),
                               ('vertical_profiles.csv',('z_m',)),
                               ('threshold_sensitivity.csv',('threshold_s-1',))]:
        # Raw records retain their actual sampling times; scalar interpolation is not a field trajectory.
        rows=[]
        for dx,p in CASES.items():
            rows.extend(dict(grid_dx_m=dx,**r) for r in read(p/filename))
        write(OUT/filename,rows)
    budget=[]
    for dx,p in CASES.items():
        budget.extend(dict(grid_dx_m=dx,**r) for r in read(p/'vorticity_budget.csv'))
    write(OUT/'vorticity_budget.csv',budget)
    domains={dx:json.loads((ROOT/('outputs/domain_extent_audit' if dx==600 else 'outputs/domain_extent_300m_audit')/'summary.json').read_text()) for dx in CASES}
    # The four calls are measured from the executed per-step control path
    # (_predictor, _project and _transport); the initialization call is excluded.
    core_tree=ast.parse((ROOT/'src/storm_dynamics/core.py').read_text(encoding='utf-8'))
    per_step_methods={'_predictor','_project','_transport'}
    calls_per_step=0
    for node in core_tree.body:
        if isinstance(node,ast.ClassDef) and node.name=='StormSimulation':
            for method in node.body:
                if isinstance(method,ast.FunctionDef) and method.name in per_step_methods:
                    calls_per_step += sum(
                        isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
                        and n.func.attr=='apply_velocity_bcs' for n in ast.walk(method)
                    )
    assert calls_per_step==4,calls_per_step
    bc_tree=ast.parse((ROOT/'src/meteorological_flow/boundary_conditions.py').read_text(encoding='utf-8'))
    damping_constants=[]
    for node in bc_tree.body:
        if isinstance(node,ast.FunctionDef) and node.name=='apply_velocity_bcs':
            for expr in ast.walk(node):
                if (isinstance(expr,ast.BinOp) and isinstance(expr.op,ast.Mult)
                        and isinstance(expr.left,ast.Constant) and expr.left.value==0.05):
                    damping_constants.append(float(expr.left.value))
    assert damping_constants==[0.05],damping_constants
    damping_fraction=damping_constants[0]
    damping={}
    for dx,seq in [(600,'diagnostic_sequence_20260905'),(300,'resolution_300m_20260909')]:
        with h5py.File(ROOT/'outputs'/seq/'sequence.h5','r') as f:
            a=f['steps'][:];weights=np.maximum(0,np.minimum(a[:,0]+a[:,1],last)-np.maximum(a[:,0],first))
            active=weights>0
            integrated_calls=float(np.sum(calls_per_step*weights[active]/a[active,1]))
            damping[dx]=dict(common_interval_s=[first,last],overlapping_steps=int(active.sum()),
                             velocity_bc_calls_per_native_step=calls_per_step,
                             maximum_fraction_per_call=damping_fraction,
                             fractional_calls=integrated_calls,
                             nominal_lowest_face_decay_rate_s_1=float(-np.log1p(-damping_fraction)*integrated_calls/(last-first)),
                             semantics='Nominal coefficient exposure only, not measured energy loss or causal influence')
    result=dict(common_times_s=common.tolist(),time_method='linear interpolation of scalar diagnostics; no extrapolation',
                analysis_geometry=expected_geometry,
                physical_metrics=stats,resolution_adequacy_metrics=adequacy,
                domain=domains,nominal_damping=damping,
                classification='Joint sensitivity to grid resolution and timestep-dependent numerical effects',
                limitations=['Different storm evolutions from analytic initialization, not paired trajectories.',
                             'Connected threshold masks share criteria but their physical sizes may differ.',
                             'Cylinder-integrated metrics use common radius and clipped 0–2 km cell volumes.',
                             'Native budget means use evolving connected masks; they are not a fixed-volume causal comparison.',
                             'Provenance initial labels start 0.326061 s apart and classify different antecedent states.',
                             'No proof that top damping is dynamically irrelevant.'])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    artifacts={}
    for path in sorted(OUT.iterdir()):
        if path.is_file() and path.name!='manifest.json':
            artifacts[path.name]={
                'bytes':path.stat().st_size,
                'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            }
    (OUT/'manifest.json').write_text(json.dumps({'artifacts':artifacts},indent=2),encoding='utf-8')


if __name__=='__main__':main()
