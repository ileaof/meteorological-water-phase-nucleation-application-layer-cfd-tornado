"""Local rain-only interventions on archived states; never advances a storm."""
from pathlib import Path
import argparse,csv,hashlib,json,sys,platform
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tests')]
from precip_microphysics import processes as proc, constants as C, thermo as th
from precip_microphysics.config import MicrophysicsConfig
from precip_microphysics.state import MicrophysicsState
from meteorological_flow.thermodynamics import theta_from_T,P0_REF
from reference_rain_evaporation_legacy import rain_evaporation as legacy


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('sequence',type=Path);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=True,parents=True)
    saved=[];dt=.5;factors=(0.,.5,.9,1.,1.1)
    with h5py.File(args.sequence,'r') as f:
        nk=int(np.searchsorted(f['grid/zc'][:],2000,side='right'));z=f['grid/zc'][:nk]
        for index,g in enumerate(f['snapshots'].values()):
            t=float(g.attrs['time_s'])
            if not 2400<=t<=3301:continue
            arrays={key:g[key][:,:,:nk] for key in ('T','rho','qv','ql','qr','qi','qs','qg','qh','p')}
            shape=arrays['T'].shape
            arrays['P']=np.broadcast_to(f['base/p0'][:nk],shape)+arrays.pop('p')
            Sw=th.saturation_ratio_water(arrays['qv'],arrays['T'],arrays['P'])
            active=np.flatnonzero((Sw<1)&(arrays['qr']>C.QSMALL))
            # Deterministic spatial coverage, with separate active and all-cell strata.
            for stratum,eligible,limit in (('active',active,512),('all',np.arange(np.prod(shape)),256)):
                if not len(eligible):continue
                ids=eligible[np.linspace(0,len(eligible)-1,min(limit,len(eligible)),dtype=int)]
                for cell in ids:
                    k=int(cell%nk)
                    saved.append(dict(index=index,time_s=t,stratum=stratum,flat_cell=int(cell),z_m=float(z[k]),
                        theta0=float(f['base/theta0'][k]),qv0=float(f['base/qv0'][k]),
                        **{name:float(a.flat[cell]) for name,a in arrays.items()}))
            print(f'local archived states t={t:.1f}',flush=True)
        steps=f['steps'][:]; selected=steps[(steps[:,0]>=2400)&(steps[:,0]<=3301)]
    write_csv(args.out/'sampled_states.csv',saved)
    a={k:np.array([r[k] for r in saved]) for k in saved[0] if k!='stratum'}
    st=MicrophysicsState(T=a['T'],P=a['P'],rho=a['rho'],qv=a['qv'],qc=a['ql'],qr=a['qr'],qi=a['qi'],qs=a['qs'],qg=a['qg'],qh=a['qh'])
    baseline=legacy(st,MicrophysicsConfig(),dt)[0].dq
    modern=proc.rain_evaporation(st,MicrophysicsConfig(),dt)[0].dq
    assert np.array_equal(modern,baseline)
    Sw=th.saturation_ratio_water(st.qv,st.T,st.P)
    lam=proc.sd.lambda_slope(st.qr,st.rho,'rain')
    rate=(2*np.pi/np.maximum(st.rho,C.TINY))*(Sw-1)/proc._diffusional_denominator(st.T,st.P,'water')*proc._ventilated_capacitance(C.N0_r,lam,proc._VA['rain'],proc._VB['rain'])
    raw=np.where(np.isfinite(rate)&(Sw<1),-rate,0.)
    theta=theta_from_T(st.T,st.P,P0_REF);loading=sum(getattr(st,k) for k in ('qc','qr','qi','qs','qg','qh'))
    tv=theta*(1+.61*st.qv-loading)
    b=9.81*((theta-a['theta0'])/a['theta0']+.61*(st.qv-a['qv0'])-loading)
    results={}; rows=[]; statistics=[]
    for factor in factors:
        tr=proc.rain_evaporation(st,MicrophysicsConfig(rain_evaporation_factor=factor),dt)
        dq=tr[0].dq if tr else np.zeros_like(baseline)
        temp=st.T-C.Lv/C.cp_d*dq;thetanew=theta_from_T(temp,st.P,P0_REF)
        tvnew=thetanew*(1+.61*(st.qv+dq)-(loading-dq))
        bnew=9.81*((thetanew-a['theta0'])/a['theta0']+.61*(st.qv+dq-a['qv0'])-(loading-dq))
        values=dict(raw_rate_kgkg_s=factor*raw,effective_rate_kgkg_s=dq/dt,dqv_kgkg=dq,dT_K=temp-st.T,
                    dtheta_K=thetanew-theta,dthetav_K=tvnew-tv,dB_ms2=bnew-b)
        results[factor]=values
        for n,r in enumerate(saved):rows.append(dict(index=r['index'],time_s=r['time_s'],stratum=r['stratum'],flat_cell=r['flat_cell'],factor=factor,dt_s=dt,**{k:float(v[n]) for k,v in values.items()}))
        mask=np.array([r['stratum']=='active' for r in saved])
        item=dict(factor=factor,active_samples=int(mask.sum()),capped_fraction=float(np.mean((dq<factor*raw*dt*(1-1e-12))[mask])) if factor else 0.)
        for key,v in values.items():
            for label,stat in zip(('p05','median','p95'),np.quantile(v[mask],[.05,.5,.95])):item[key+'_'+label]=float(stat)
        statistics.append(item)
    write_csv(args.out/'local_responses.csv',rows);write_csv(args.out/'local_statistics.csv',statistics)
    fig,axs=plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
    keys=('effective_rate_kgkg_s','dqv_kgkg','dT_K','dtheta_K','dthetav_K','dB_ms2')
    for ax,key in zip(axs.flat,keys):
        med=[r[key+'_median'] for r in statistics];lo=[r[key+'_p05'] for r in statistics];hi=[r[key+'_p95'] for r in statistics]
        ax.plot(factors,med,'o-');ax.fill_between(factors,lo,hi,alpha=.2);ax.set(xlabel='Fator de evaporação da chuva',ylabel=key);ax.grid(alpha=.2)
    fig.suptitle('Intervenção local de 0,5 s: mediana e P05–P95 dos estados ativos; sem evolução da tempestade')
    fig.savefig(args.out/'local_sensitivity.png');plt.close(fig)
    summary=dict(no_storm_run=True,dt_local_s=dt,factors_tested=factors,perturbation_for_future_experiment_selected=False,
        sampled_states=len(saved),snapshots=len(set(r['index'] for r in saved)),method='512 active + 256 all-cell evenly spaced flat indices per snapshot, deterministic; overlaps possible, not a population-weighted sample',
        archived_dt_range_s=[float(selected[:,1].min()),float(selected[:,1].max())],
        baseline_identity={'max_error':float(np.max(abs(modern-baseline))),'rms_error':float(np.sqrt(np.mean((modern-baseline)**2))),'bitwise_equal':True},
        state_ranges={k:[float(a[k].min()),float(a[k].max())] for k in ('T','P','rho','qv','qr','z_m')},statistics=statistics,
        relative_to_f1={str(factor):{key:float(np.median((results[factor][key]-results[1.][key])[mask])) for key in keys} for factor in (.9,1.1)},
        python=sys.version,platform=platform.platform(),source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),ROOT/'src/precip_microphysics/processes.py',ROOT/'src/precip_microphysics/config.py']})
    (args.out/'local_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k in ('sampled_states','baseline_identity','relative_to_f1','archived_dt_range_s')},indent=2))


if __name__=='__main__':main()
