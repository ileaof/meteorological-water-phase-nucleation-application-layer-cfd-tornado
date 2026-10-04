"""Authorized 600 m replay with a mandatory short gate before continuation."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from run_causal_evaporation_branch import PROGNOSTIC, restore
from run_diagnostic_sequence import build
from storm_dynamics.muscl_capture import DirectMusclCapture, TERMS

ENDPOINTS = (0,1,3,4,8,9,11,14,15,17)
RADII = (1200.,2400.,4200.,6000.)
MODES = ('moving_end','moving_symmetric','fixed_initial')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str), encoding='utf-8')


def masks_for_block(x, y, track, a, b):
    def disc(index, radius):
        row = track[index]
        return (((x[:,None]-row['center_x_m'])**2+(y[None,:]-row['center_y_m'])**2)<=radius**2).astype(float)
    masks, pairs, descriptions = [], [], []
    for radius in RADII:
        start, end, fixed = disc(a,radius), disc(b,radius), disc(0,radius)
        for mode in MODES:
            ma,mb = (fixed,fixed) if mode=='fixed_initial' else (start,end)
            symmetric = mode=='moving_symmetric'
            masks.append((ma+mb)*.5 if symmetric else mb)
            pairs.append((ma,mb,symmetric))
            descriptions.append(dict(radius_m=radius,mask_mode=mode))
    return np.stack(masks),pairs,descriptions


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=ROOT/'outputs/muscl_direct_capture_20260920')
    parser.add_argument('--device',choices=('cpu','gpu'),default='gpu')
    parser.add_argument('--pilot-only',action='store_true')
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=False)
    sequence=ROOT/'outputs/diagnostic_sequence_20260905/sequence.h5'
    provenance=ROOT/'outputs/vorticity_provenance_long_v4_2790_3300/provenance_long_v4.h5'
    track_path=ROOT/'outputs/resolution_600m_comparison_audit_20260909/vortex_timeseries.csv'
    prior_path=ROOT/'outputs/les_discrete_balance_20260915_v2/terms_by_block.csv'
    with track_path.open(encoding='utf-8') as stream:
        track=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(stream)]
    seq=h5py.File(sequence,'r'); prov=h5py.File(provenance,'r')
    first=seq['snapshots/00093']; last=seq['snapshots/00110']
    start_step=int(first.attrs['step']); end_step=int(last.attrs['step'])
    start_t=float(first.attrs['time_s']); end_t=float(last.attrs['time_s'])
    all_steps=seq['steps'][:]
    schedule=all_steps[(all_steps[:,2]>=start_step)&(all_steps[:,2]<end_step)]
    assert len(schedule)==1015
    assert np.array_equal(schedule[:,2],np.arange(start_step,end_step))
    assert abs(schedule[0,0]-start_t)<2e-8
    assert np.max(abs(schedule[:-1,0]+schedule[:-1,1]-schedule[1:,0]))<2e-8
    assert abs(schedule[-1,0]+schedule[-1,1]-end_t)<2e-8
    x,y,z=[seq['grid/'+n][:] for n in ('xc','yc','zc')]
    block_steps=[int(prov[f'snapshots/{i:05d}'].attrs['step']) for i in ENDPOINTS]
    snapshots={int(g.attrs['step']):g.name for g in seq['snapshots'].values()
               if start_step<=int(g.attrs['step'])<=end_step}
    estimate=len(schedule)*len(TERMS)*12*17*3*8+2*1024**2
    free=shutil.disk_usage(args.out).free
    if free<1024**3+estimate:
        raise RuntimeError('Insufficient storage reserve')
    metadata=dict(status='preflight',sequence=str(sequence),provenance=str(provenance),
        command=sys.argv,grid=[120,120,48],start_time_s=start_t,end_time_target_s=end_t,
        steps_target=1015,pilot_steps=block_steps[1]-start_step,
        estimated_output_bytes=estimate,disk_free_bytes_before=free,
        schedule_sha256=hashlib.sha256(np.asarray(schedule,dtype='<f8').tobytes()).hexdigest(),
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for directory in ('src','scripts')
                       for p in (ROOT/directory).rglob('*.py')},
        input_tables_sha256={str(p.relative_to(ROOT)):sha(p) for p in (track_path,prior_path)},
        flux_semantics='native routines reevaluated on their unmodified stage input; checked against actual updates',
        units=dict(flux='m2/s2',velocity_increment='m/s',curl_increment='1/s',circulation_increment='m2/s'),
        tolerance=dict(relative=1e-10,velocity_rms_floor_m_s=1e-12,
                       curl_rms_floor_s_1=1e-12,circulation_rms_floor_m2_s=1.,time_absolute_s=2e-8),
        cpu_numpy_version=np.__version__,python=sys.version)
    write_json(args.out/'metadata.json',metadata)
    print(f'Estimated raw output {estimate/1024**2:.1f} MiB; first {metadata["pilot_steps"]} steps are mandatory pilot.',flush=True)
    started=time.perf_counter()
    plain,_=build(args.device,120,48,end_t,1.)
    observed,_=build(args.device,120,48,end_t,1.)
    restore(plain,sequence,93);restore(observed,sequence,93)
    xp=observed.grid.xp;to=observed.grid.backend.to_cpu
    masks,pairs,descriptions=masks_for_block(x,y,track,ENDPOINTS[0],ENDPOINTS[1])
    capture=DirectMusclCapture(observed,masks)
    observed.diagnostic_observer=capture
    original_checks=[];old_label_checks=[];comparison=[]
    source_digest=hashlib.sha256()
    consumed=[]
    def reference(group,name,selection=...):
        arr=np.asarray(group[name][selection])
        tag=dict(path=group[name].name,selection=repr(selection),shape=arr.shape,dtype=str(arr.dtype))
        source_digest.update(json.dumps(tag,sort_keys=True).encode())
        source_digest.update(np.ascontiguousarray(arr).tobytes());consumed.append(tag)
        return arr
    def check_original(sim,step):
        if step not in snapshots:return
        group=seq[snapshots[step]]
        row=dict(step=step,time_s=sim.t,fields={})
        for name in PROGNOSTIC:
            current=to(getattr(sim.state,name)); original=reference(group,name)
            equal=bool(np.array_equal(current,original))
            row['fields'][name]=dict(bitwise=equal,max_abs=float(np.max(abs(current-original))))
        original_checks.append(row)
        if not all(v['bitwise'] for v in row['fields'].values()):
            raise RuntimeError('Original snapshot mismatch: '+json.dumps(row))
    check_original(observed,start_step)
    file=h5py.File(args.out/'capture.h5','x')
    file.attrs['schema']='storm-muscl-direct-capture-v1'
    file.attrs['status']='running';file.attrs['terms']=json.dumps(TERMS)
    file.attrs['masks']=json.dumps(descriptions)
    file.attrs['stat_columns']='signed_integral,spatial_abs_integral,adjoint_integral'
    file.attrs['balance_columns']='inventory_start,inventory_end,mask_motion'
    file.attrs['balance_quantities']='total,les_label'
    file.create_dataset('schedule',data=schedule)
    file.create_dataset('levels_m',data=z[:17])
    data=file.create_dataset('step_statistics',shape=(1015,len(TERMS),12,17,3),dtype='f8',
        chunks=(1,len(TERMS),12,17,3),compression='gzip',compression_opts=1,fillvalue=np.nan)
    blockdata=file.create_dataset('block_statistics',shape=(9,len(TERMS),12,17,3),dtype='f8',fillvalue=np.nan)
    balancedata=file.create_dataset('block_balances',shape=(9,12,2,17,3),dtype='f8',fillvalue=np.nan)
    blockgroup=file.create_group('blocks')
    with prior_path.open(encoding='utf-8') as stream:
        oldrows=[r for r in csv.DictReader(stream) if int(r['grid_dx_m'])==600]
    oldindex={(int(r['block']),float(r['radius_m']),r['mask_mode'],float(r['level_m']),r['quantity'],r['term']):float(r['signed_m2_s']) for r in oldrows}
    block=0;completed_steps=0;failed=None;pilot_pass=False
    try:
        for number,(t0,dt,step) in enumerate(schedule,start=1):
            for sim in (plain,observed):
                if sim.step!=int(step) or abs(sim.t-t0)>2e-8:raise RuntimeError('Schedule mismatch before step')
                sim._step(float(dt));sim.step+=1;sim.t=float(sim.state.t)
                if abs(sim.t-(t0+dt))>2e-8:raise RuntimeError('Schedule mismatch after step')
            for name in PROGNOSTIC:
                if not bool(xp.array_equal(getattr(plain.state,name),getattr(observed.state,name))):
                    raise RuntimeError('Observer neutrality failed: '+name)
            check_original(observed,observed.step)
            data[number-1]=capture.last_record
            completed_steps=number
            file.attrs['completed_steps']=number
            if observed.step==block_steps[block+1]:
                stats,balances=capture.finish_block(pairs)
                blockdata[block]=stats;balancedata[block]=balances
                group=blockgroup.create_group(f'{block:02d}')
                group.attrs['start_step']=block_steps[block]
                group.attrs['end_step']=observed.step
                group.create_dataset('masks',data=masks,compression='gzip')
                pgroup=prov[f'snapshots/{ENDPOINTS[block+1]:05d}']
                for quantity,name in enumerate(('omega_total','omega_source/les')):
                    old=reference(pgroup,name,(2,slice(None),slice(None),slice(0,17)))
                    current=capture.q()[quantity]
                    capture.check('old_provenance_'+str(quantity),current,xp.asarray(old),1e-12,'1/s')
                    old_label_checks.append(dict(block=block,quantity=name,bitwise=bool(np.array_equal(to(current),old))))
                for m,desc in enumerate(descriptions):
                    for term,quantity,oldterm in (('total_les','total','les'),('total_muscl','total','advection'),
                        ('les_injection','les_label','local_les'),('les_muscl','les_label','inferred_muscl_label')):
                        expected=np.array([oldindex[(block,desc['radius_m'],desc['mask_mode'],float(level),quantity,oldterm)] for level in z[:17]])
                        measured=stats[TERMS.index(term),m,:,0]
                        capture.check('previous_balance_'+term,xp.asarray(measured),xp.asarray(expected),1.,'m2/s')
                        comparison.append(dict(block=block,**desc,term=term,max_abs_m2_s=float(np.max(abs(measured-expected)))))
                file.attrs['completed_blocks']=block+1;file.flush()
                if block==0:
                    pilot_pass=True
                    write_json(args.out/'pilot_gate.json',dict(status='PASS',steps=number,
                        original_snapshots_checked=len(original_checks),gates=capture.gates,
                        wall_clock_s=time.perf_counter()-started,
                        estimated_full_wall_s=(time.perf_counter()-started)*1015/number,
                        neutral_bitwise_all_12=True))
                    print('PILOT PASS: original snapshots, bitwise neutrality and all reconstruction gates passed.',flush=True)
                    if args.pilot_only:break
                block+=1
                if block<9:
                    masks,pairs,_=masks_for_block(x,y,track,ENDPOINTS[block],ENDPOINTS[block+1])
                    capture.set_masks(masks)
            if number==1 or number%10==0 or number==1015:
                elapsed=time.perf_counter()-started
                print(f'step {number}/1015 t={observed.t:.6f} block={block} wall={elapsed:.1f}s pilot={pilot_pass}',flush=True)
    except Exception as exc:
        failed=f'{type(exc).__name__}: {exc}'
        raise
    finally:
        status='failed' if failed else ('pilot_complete' if args.pilot_only else 'complete')
        file.attrs['status']=status;file.flush();file.close()
        seq.close();prov.close()
        metadata.update(status=status,failure=failed,completed_steps=completed_steps,
            final_time_s=observed.t,pilot_pass=pilot_pass,wall_clock_s=time.perf_counter()-started,
            gates=capture.gates,original_snapshots=original_checks,old_label_checks=old_label_checks,
            prior_balance_comparison=comparison,neutrality_bitwise_all_12=(failed is None),
            input_arrays_sha256=source_digest.hexdigest(),input_selections=consumed)
        write_json(args.out/'metadata.json',metadata)
        digest=hashlib.sha256()
        with (args.out/'capture.h5').open('rb') as stream:
            while chunk:=stream.read(8*1024**2):digest.update(chunk)
        write_json(args.out/'manifest.json',dict(capture_sha256=digest.hexdigest(),
            capture_bytes=(args.out/'capture.h5').stat().st_size,
            metadata_sha256=sha(args.out/'metadata.json'),source_sha256=metadata['source_sha256']))
        print(f'{status}: {completed_steps} steps, wall={metadata["wall_clock_s"]:.1f}s',flush=True)


if __name__=='__main__':main()
