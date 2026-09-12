"""Reproduce the existing idealized stress-drag parent, with passive capture.

No historical cache is imported. Physics settings are copied from
scratchpad/tornado_intensity_L_gpu.py's parent (no nests). Adaptive timesteps
are not adjusted to hit output times; actual synchronized times are saved.
"""
from pathlib import Path
import argparse
import dataclasses
import hashlib
import json
import platform
import shutil
import subprocess
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.soundings import bunkers_storm_motion
from meteorological_flow.base_state import BaseState
from storm_dynamics.diagnostic_capture import DiagnosticCapture


def build(device='gpu',nx=120,nz=48,duration=3900.,rain_evaporation_factor=1.0,
          *,ny=None,Lx_m=72000.,Ly_m=None,Lz_m=15000.,high_top_m=None,
          damping_faces=None,reference_T=None,reference_qv=None):
    ny=nx if ny is None else ny
    Ly_m=Lx_m if Ly_m is None else Ly_m
    cfg=build_storm_config(preset='storm',nx=nx,ny=ny,nz=nz,Lx=Lx_m,Ly=Ly_m,Lz=Lz_m,
                          duration=duration,dt_max=3.,drag=True,z_stretch=1.05,C_s=.20,
                          hodograph_kind='quarter_circle',U_max=30.,device=device)
    cfg.sim.physics.bubble_dtheta=5.
    cfg.sim.physics.rain_evaporation_factor=float(rain_evaporation_factor)
    if reference_T is not None:
        cfg.sim.physics.T_ref=float(reference_T)
    if damping_faces is not None:
        cfg.sim.boundaries.damping_faces=int(damping_faces)
    cfg.dyn.drag.stress_divergence=True; cfg.dyn.drag.surface_layer_depth_m=150.
    cfg.dyn.drag.use_log_law=True; cfg.dyn.drag.roughness_length_m=.1
    if high_top_m is not None:
        from storm_dynamics.top_boundary import configure_high_top
        cfg=configure_high_top(cfg,float(high_top_m))
    initial=StormSimulation(cfg); b=initial.base; cx,cy=bunkers_storm_motion(b)
    base=BaseState(zc=b.zc,theta0=b.theta0,qv0=b.qv0,p0=b.p0,T0=b.T0,rho0=b.rho0,
                   u0=b.u0-cx,v0=b.v0-cy)
    del initial
    sim=StormSimulation(cfg,base=base)
    if reference_qv is not None:
        sim.qv_ref=float(reference_qv)
    return sim,(cx,cy)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--duration',type=float,default=3900.)
    parser.add_argument('--interval',type=float,default=30.)
    parser.add_argument('--device',choices=['cpu','gpu'],default='gpu')
    parser.add_argument('--nx',type=int,default=120); parser.add_argument('--ny',type=int,default=None)
    parser.add_argument('--nz',type=int,default=48)
    parser.add_argument('--Lx-m',type=float,default=72000.); parser.add_argument('--Ly-m',type=float,default=None)
    parser.add_argument('--Lz-m',type=float,default=15000.)
    parser.add_argument('--rain-evaporation-factor',type=float,default=1.0)
    parser.add_argument('--capture-start',type=float,default=0.,help='First native time to begin capture')
    args=parser.parse_args()
    if not 0 <= args.capture_start < args.duration:
        parser.error('capture-start must be nonnegative and less than duration')
    args.out.mkdir(parents=True,exist_ok=False)
    start=time.time(); sim,motion=build(
        device=args.device,nx=args.nx,ny=args.ny,nz=args.nz,
        Lx_m=args.Lx_m,Ly_m=args.Ly_m,Lz_m=args.Lz_m,
        duration=args.duration,rain_evaporation_factor=args.rain_evaporation_factor,
    )
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for directory in ('src','scripts') for p in (ROOT/directory).rglob('*.py')}
    metadata=dict(config=dataclasses.asdict(sim.scfg),source_sha256=hashes,
                  git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  git_diff=subprocess.check_output(['git','diff'],cwd=ROOT,text=True),
                  command=sys.argv,python=sys.version,numpy=np.__version__,platform=platform.platform(),
                  backend=sim.backend.device_info(),frame='storm-relative, fixed grid',
                  storm_motion_ground_ms=motion,initial_condition='fresh analytic sounding and warm bubble; no cache',
                  physics_source='scratchpad/tornado_intensity_L_gpu.py parent configuration',
                  output_interval_target_s=args.interval,prognostic_step_clipping='only at final integration time',
                  budget_top_m=2000.,capture_start_target_s=args.capture_start,status='running')
    (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str),encoding='utf-8')
    capture=None
    if args.capture_start == 0:
        capture=DiagnosticCapture(sim,args.out/'sequence.h5',metadata,args.interval)
        sim.diagnostic_observer=capture
    last_report=-1.; completed=False
    try:
        while sim.t<args.duration-1e-9:
            dt=min(float(sim._dt()),args.duration-sim.t)
            sim._step(dt); sim.step+=1; sim.t=float(sim.state.t)
            if capture is None and sim.t >= args.capture_start:
                capture=DiagnosticCapture(sim,args.out/'sequence.h5',metadata,args.interval)
                sim.diagnostic_observer=capture
            if sim.t-last_report>=30. or sim.t>=args.duration:
                if not bool(sim.grid.xp.isfinite(sim.state.w).all()): raise RuntimeError('Nonfinite velocity')
                if shutil.disk_usage(args.out).free<3*1024**3: raise RuntimeError('Capture stopped: less than 3 GiB free')
                message=f't={sim.t:.3f}s step={sim.step} frames={capture.count if capture else 0} wall={time.time()-start:.1f}s wmax={float(sim.state.w.max()):.3f}'
                print(message,flush=True)
                with (args.out/'progress.log').open('a',encoding='utf-8') as stream: stream.write(message+'\n')
                last_report=sim.t
        completed=True
    finally:
        if capture is not None:
            capture.close(sim,completed)
        metadata.update(status='complete' if completed else 'interrupted',final_time_s=sim.t,
                        steps=sim.step,wall_clock_s=time.time()-start,frames=capture.count if capture else 0)
        (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str),encoding='utf-8')


if __name__=='__main__': main()
