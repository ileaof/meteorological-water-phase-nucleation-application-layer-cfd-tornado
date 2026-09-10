"""Restart-controlled rain-evaporation branch using the archived native dt schedule.

This is a pilot intervention beginning at one archived synchronized state.  It
does not claim to reproduce an intervention active since model initialization.
"""
from pathlib import Path
import argparse,csv,dataclasses,hashlib,json,platform,subprocess,sys,time
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
from run_diagnostic_sequence import build
from storm_dynamics.diagnostic_capture import DiagnosticCapture,kinematics


PROGNOSTIC=('u','v','w','p','theta','qv','ql','qi','qr','qs','qg','qh')


def restore(sim,sequence,snapshot):
    with h5py.File(sequence,'r') as f:
        group=f[f'snapshots/{snapshot:05d}']; xp=sim.grid.xp
        for name in PROGNOSTIC:
            setattr(sim.state,name,xp.asarray(group[name][:]))
        sim.state.p_dyn=xp.asarray(group['p_dyn'][:])
        sim.state.t=float(group.attrs['time_s']);sim.t=sim.state.t;sim.step=int(group.attrs['step'])
        sim.state.surface_precip={name:xp.zeros((sim.grid.nx,sim.grid.ny)) for name in ('rain','snow','graupel','hail')}
        sim.state.diagnose(sim.cfg)
        for name in ('xc','yc','zc','xf','yf','zf'):
            if not np.array_equal(sim.grid.backend.to_cpu(getattr(sim.grid,name)),f['grid/'+name][:]):
                raise RuntimeError('Restart grid mismatch: '+name)
        for name in ('theta0','qv0','p0','T0','rho0','u0','v0'):
            if not np.array_equal(np.asarray(getattr(sim.base,name)),f['base/'+name][:]):
                raise RuntimeError('Restart base-state mismatch: '+name)
        return sim.t,sim.step


class Scalars:
    def __init__(self,sim):
        g=sim.grid;xp=g.xp
        dz=g.dz_c if getattr(g,'stretched',False) else xp.full(g.nz,g.dz)
        self.weight=sim.rho0_c[None,None,:]*float(g.dx)*float(g.dy)*dz[None,None,:]
        self.mass={};self.energy={};self.sim=sim
        sim.coupler.scheme.process_observer=self.observe

    def observe(self,name,dq,latent_per_dq):
        value=(dq*self.weight).sum()
        self.mass[name]=self.mass.get(name,0.0)+value
        self.energy[name]=self.energy.get(name,0.0)+value*float(latent_per_dq)*1005.0

    def row(self):
        sim=self.sim;g=sim.grid;st=sim.state;xp=g.xp;to=g.backend.to_cpu
        loading=st.ql+st.qi+st.qr+st.qs+st.qg+st.qh
        tv=st.theta*(1+.61*st.qv-loading);tv0=sim.theta0_field*(1+.61*sim.qv0_field)
        anomaly=tv-tv0;B=sim.cfg.flow.gravity*((st.theta-sim.theta0_field)/sim.theta0_field+.61*(st.qv-sim.qv0_field)-loading)
        k=int(np.argmin(abs(to(g.zc)-100))); cold=anomaly[:,:,k]<-1.; area=float(g.dx*g.dy)/1e6
        z=to(g.zc); depth=xp.where(anomaly<-1.,xp.asarray(z)[None,None,:],0.).max(axis=2)
        # Periodic horizontal derivative, used only for scalar screening metrics.
        bx=(xp.roll(B,-1,axis=0)-xp.roll(B,1,axis=0))/(2*float(g.dx))
        by=(xp.roll(B,-1,axis=1)-xp.roll(B,1,axis=1))/(2*float(g.dy));bgrad=xp.hypot(bx,by)
        inds=np.where(z<=500)[0]; kin=kinematics(tuple(to(getattr(st,n)[:,:,:len(inds)+(n=='w')]) for n in ('u','v','w')),
            (to(g.xc),to(g.yc),z[inds]))
        zeta=kin['omega'][2];ii=np.unravel_index(np.nanargmax(zeta),zeta.shape)
        xi,eta=kin['omega'][0][ii],kin['omega'][1][ii]
        wc=(to(st.w[:,:,:len(inds)])+to(st.w[:,:,1:len(inds)+1]))/2
        wx,wy,_=np.gradient(wc,to(g.xc),to(g.yc),z[inds],edge_order=2)
        gx,gy=wx[ii],wy[ii];norm=np.hypot(xi,eta)*np.hypot(gx,gy)
        result=dict(time_s=float(sim.t),step=sim.step,cold_area_m1_km2=float(to(cold.sum()))*area,
            cold_area_m05_km2=float(to((anomaly[:,:,k]<-.5).sum()))*area,cold_area_m2_km2=float(to((anomaly[:,:,k]<-2).sum()))*area,
            theta_v_min_K=float(to(anomaly.min())),theta_v_cold_mean_K=float(to(anomaly[:,:,k][cold].mean())) if bool(to(cold.any())) else np.nan,
            B_min_ms2=float(to(B.min())),B_cold_mean_ms2=float(to(B[:,:,k][cold].mean())) if bool(to(cold.any())) else np.nan,
            cold_depth_mean_m=float(to(depth[cold].mean())) if bool(to(cold.any())) else 0.,cold_depth_max_m=float(to(depth.max())),
            B_gradient_max_s2=float(to(bgrad.max())),B_gradient_cold_mean_s2=float(to(bgrad[:,:,k][cold].mean())) if bool(to(cold.any())) else np.nan,
            w_max_ms=float(to(st.w.max())),total_condensate_kg=float(to((loading*self.weight).sum())),
            zeta_low_max_s=float(zeta[ii]),zeta_peak_x_m=float(to(g.xc)[ii[0]]),zeta_peak_y_m=float(to(g.yc)[ii[1]]),zeta_peak_z_m=float(z[ii[2]]),
            xi_s=float(xi),eta_s=float(eta),wx_s=float(gx),wy_s=float(gy),tilting_s2=float(xi*gx+eta*gy),stretching_s2=float(kin['stretching'][ii]),
            omega_orientation_deg=float(np.degrees(np.arctan2(eta,xi))),gradw_orientation_deg=float(np.degrees(np.arctan2(gy,gx))),
            omega_gradw_angle_deg=float(np.degrees(np.arccos(np.clip((xi*gx+eta*gy)/norm,-1,1)))) if norm else np.nan)
        for name,value in self.mass.items():result['cumulative_'+name+'_mass_kg']=float(to(value))
        for name,value in self.energy.items():result['cumulative_'+name+'_latent_J']=float(to(value))
        for name,value in st.surface_precip.items():result['surface_'+name+'_kg']=float(to(value.sum()))*float(g.dx*g.dy)
        return result


def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sequence',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--factor',type=float,required=True);p.add_argument('--start-snapshot',type=int,default=79);p.add_argument('--end-snapshot',type=int,required=True)
    p.add_argument('--device',choices=('cpu','gpu'),default='gpu');p.add_argument('--capture',action='store_true');p.add_argument('--validate-control',action='store_true')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False);start_wall=time.time()
    with h5py.File(a.sequence,'r') as f:
        end_time=float(f[f'snapshots/{a.end_snapshot:05d}'].attrs['time_s']);saved={round(float(g.attrs['time_s']),9):g.name for g in f['snapshots'].values()}
        steps=f['steps'][:];historical_metadata=json.loads(f.attrs['metadata'])
    sim,motion=build(a.device,120,48,end_time,a.factor);start_time,start_step=restore(sim,a.sequence,a.start_snapshot)
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for d in ('src','scripts') for p in (ROOT/d).rglob('*.py')}
    metadata=dict(case=a.out.name,rain_evaporation_factor=a.factor,intervention_start_s=start_time,end_time_s=end_time,restart_snapshot=a.start_snapshot,
        restart_semantics='all prognostic arrays restored; surface precipitation reset for post-intervention accumulation only',
        timestep_semantics='exact archived CONTROL dt schedule imposed on every branch',config=dataclasses.asdict(sim.scfg),config_sha256='',source_sha256=hashes,
        git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),git_diff=subprocess.check_output(['git','diff'],cwd=ROOT,text=True),
        command=sys.argv,python=sys.version,numpy=np.__version__,platform=platform.platform(),backend=sim.backend.device_info(),seed=sim.cfg.random_seed,
        storm_motion_ground_ms=motion,historical_sequence=str(a.sequence.resolve()),historical_config=historical_metadata['config'])
    metadata['config_sha256']=hashlib.sha256(json.dumps(metadata['config'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    (a.out/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str),encoding='utf-8')
    scalar=Scalars(sim);capture=DiagnosticCapture(sim,a.out/'sequence.h5',metadata,30.,full_column=False) if a.capture else None
    # DiagnosticCapture installs its own compatible process observer when active.
    if capture:
        scalar=None
        # Preserve the historical absolute 30-s output clock rather than
        # re-anchoring it at the slightly-late restart snapshot.
        capture.next_output=float(np.ceil(start_time/30.)*30.)
        sim.diagnostic_observer=capture
    rows=[];validation=[];next_record=start_time
    schedule=steps[(steps[:,2]>=start_step)&(steps[:,0]<end_time-1e-8)]
    for t0,dt,step_index in schedule:
        if abs(sim.t-t0)>2e-8:raise RuntimeError(f'dt schedule mismatch {sim.t} vs {t0}')
        own_dt=float(sim._dt());sim._step(float(dt));sim.step+=1;sim.t=float(sim.state.t)
        if float(dt)>own_dt*(1+1e-10):
            with (a.out/'cfl_exceedances.log').open('a') as f:f.write(f'{t0},{dt},{own_dt}\n')
        if sim.t>=next_record-1e-8:
            if scalar:rows.append(scalar.row())
            next_record+=30.
        key=round(sim.t,9)
        if a.validate_control and key in saved:
            with h5py.File(a.sequence,'r') as f:
                g=f[saved[key]];entry={'time_s':sim.t,'bitwise':True,'max_error':0.}
                for name in PROGNOSTIC:
                    current=sim.grid.backend.to_cpu(getattr(sim.state,name));reference=g[name][:]
                    entry['bitwise']&=bool(np.array_equal(current,reference));entry['max_error']=max(entry['max_error'],float(np.max(abs(current-reference))))
                validation.append(entry)
        if sim.t>=end_time-1e-8:break
    if scalar:rows.append(scalar.row());write_csv(a.out/'screening_metrics.csv',rows)
    if capture:capture.close(sim,completed=abs(sim.t-end_time)<1e-7)
    metadata.update(final_time_s=sim.t,steps_executed=sim.step-start_step,wall_clock_s=time.time()-start_wall,control_validation=validation,
                    archived_dt_min_s=float(schedule[:,1].min()),archived_dt_max_s=float(schedule[:,1].max()),status='complete')
    (a.out/'metadata.json').write_text(json.dumps(metadata,indent=2,default=str),encoding='utf-8')
    print(json.dumps({k:metadata[k] for k in ('case','rain_evaporation_factor','final_time_s','steps_executed','wall_clock_s','control_validation')},indent=2),flush=True)


if __name__=='__main__':main()
