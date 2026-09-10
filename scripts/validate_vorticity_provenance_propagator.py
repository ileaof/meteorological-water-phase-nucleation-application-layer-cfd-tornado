"""Reduced long-memory and directional validation of the provenance propagator."""
from __future__ import annotations
import json, hashlib, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from meteorological_flow.advection import _minmod
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.momentum import momentum_advection_tendency
from storm_dynamics.vorticity_provenance import (_frozen_slope,curl_native,
    frozen_momentum_advection_tendency,frozen_muscl_advection_tendency)


def sim(nx=32,ny=8,nz=6):
    return StormSimulation(build_storm_config(nx=nx,ny=ny,nz=nz,Lx=32000.,Ly=8000.,Lz=4000.,z_stretch=1.,duration=1.,device="cpu"))


def rms(a): return float(np.sqrt(np.mean(np.asarray(a)**2)))


def one_dimensional():
    s=sim(64,4,4); g=s.grid; x=np.arange(g.nx,dtype=float)
    patterns={"smooth":np.sin(2*np.pi*x/g.nx),"pulse":((x>=18)&(x<34)).astype(float),
              "gradient":(x-g.nx/2)/g.nx,"sign":np.where(x<g.nx/2,-1.,1.),
              "limited":np.sin(2*np.pi*x/g.nx)+.35*np.sin(10*np.pi*x/g.nx)}
    vel=(np.full(g.u_shape,12.),np.zeros(g.v_shape),np.zeros(g.w_shape)); dt=.35*g.dx/12.; out={}
    for name,q0 in patterns.items():
        q=np.broadcast_to(q0[:,None,None],(g.nx,g.ny,g.nz)).copy(); a=.55*q+.08*np.roll(q,3,0); b=q-a
        initial_max=float(np.max(abs(q))); max_closure=0.; kappas=[]
        for n in range(1600):
            tq=frozen_muscl_advection_tendency(q,q,vel,g); ta=frozen_muscl_advection_tendency(a,q,vel,g); tb=frozen_muscl_advection_tendency(b,q,vel,g)
            q+=dt*tq; a+=dt*ta; b+=dt*tb
            if n%40==0:
                max_closure=max(max_closure,float(np.max(abs(q-a-b))))
                kappas.append((rms(a)+rms(b))/max(rms(q),1e-30))
        out[name]={"crossing_times":1600*dt*12/(g.nx*g.dx),"max_abs_closure":max_closure,
                   "max_amplitude_ratio":float(np.max(abs(q)))/max(initial_max,1e-30),
                   "kappa_initial":kappas[0],"kappa_final":kappas[-1],"kappa_max":max(kappas)}
    # FTCS reference on the same Fourier mode.
    q=np.sin(2*np.pi*9*x/g.nx); q0=np.linalg.norm(q); c=.45
    for _ in range(300): q-=.5*c*(np.roll(q,-1)-np.roll(q,1))
    out["old_centered_euler"]={"steps":300,"amplification":float(np.linalg.norm(q)/q0)}
    return out


def epsilon_study():
    s=sim(); g=s.grid; x=np.arange(g.nx,dtype=float)[:,None,None]
    total=(x+.02*x*x)*np.ones((1,g.ny,g.nz)); d=np.sin(.17*x)*np.ones((1,g.ny,g.nz))
    exact=_frozen_slope(d,total,0,g,True); base=_minmod(total-np.roll(total,1,0),np.roll(total,-1,0)-total); rows=[]
    for eps in (1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9):
        p=total+eps*d; fd=(_minmod(p-np.roll(p,1,0),np.roll(p,-1,0)-p)-base)/eps
        rows.append({"epsilon":eps,"max_abs_error":float(np.max(abs(fd-exact))),"rms_error":rms(fd-exact)})
    return rows


def temporal_refinement():
    s=sim(64,4,4); g=s.grid; x=np.arange(g.nx,dtype=float)
    q0=np.broadcast_to(np.sin(2*np.pi*x/g.nx)[:,None,None],(g.nx,g.ny,g.nz)).copy()
    vel=(np.full(g.u_shape,12.),np.zeros(g.v_shape),np.zeros(g.w_shape))
    def evolve(dt,steps):
        q=q0.copy(); a=.55*q+.08*np.roll(q,3,0); b=q-a
        for _ in range(steps):
            tq=frozen_muscl_advection_tendency(q,q,vel,g); ta=frozen_muscl_advection_tendency(a,q,vel,g); tb=frozen_muscl_advection_tendency(b,q,vel,g)
            q+=dt*tq; a+=dt*ta; b+=dt*tb
        return q,a,b
    dt=.10*g.dx/12.; coarse=evolve(dt,700); fine=evolve(dt/2,1400)
    return {"same_final_time_s":700*dt,"coarse_dt_s":dt,"fine_dt_s":dt/2,
            "total_relative_rms":rms(coarse[0]-fine[0])/rms(fine[0]),
            "label_1_relative_rms":rms(coarse[1]-fine[1])/rms(fine[1]),
            "label_2_relative_rms":rms(coarse[2]-fine[2])/rms(fine[2])}


def long_cgrid():
    s=sim(10,10,6); g=s.grid; rng=np.random.default_rng(9)
    total=(s.state.u.copy(),s.state.v.copy(),s.state.w.copy())
    a=tuple(.5*q+.001*rng.normal(size=q.shape) for q in total); b=tuple(q-a[i] for i,q in enumerate(total)); dt=.005
    rows=[]; max_close=0.
    for n in range(2500):
        tt=frozen_momentum_advection_tendency(total,total,g); ta=frozen_momentum_advection_tendency(a,total,g); tb=frozen_momentum_advection_tendency(b,total,g)
        total=tuple(q+dt*d for q,d in zip(total,tt)); a=tuple(q+dt*d for q,d in zip(a,ta)); b=tuple(q+dt*d for q,d in zip(b,tb))
        if n%50==0 or n==2499:
            ot=curl_native(total,g); oa=curl_native(a,g); ob=curl_native(b,g); res=ot-oa-ob; max_close=max(max_close,float(np.max(abs(res))))
            k=(rms(oa)+rms(ob))/max(rms(ot),1e-30); rows.append({"step":n+1,"kappa_omega":k,"closure_max":float(np.max(abs(res)))})
    return {"steps":2500,"dt":dt,"max_closure":max_close,"kappa_initial":rows[0]["kappa_omega"],
            "kappa_final":rows[-1]["kappa_omega"],"kappa_max":max(r["kappa_omega"] for r in rows),"history":rows}


def full_jacobian_distinction():
    s=sim(12,10,6); g=s.grid
    s.state.u += .8*np.sin(2*np.pi*g.xf[:,None,None]/(g.nx*g.dx))*np.cos(2*np.pi*g.yc[None,:,None]/(g.ny*g.dy))
    s.state.v += -.8*(g.ny*g.dy/(g.nx*g.dx))*np.cos(2*np.pi*g.xc[:,None,None]/(g.nx*g.dx))*np.sin(2*np.pi*g.yf[None,:,None]/(g.ny*g.dy))
    s.state.w[:,:,1:-1] += .1*np.sin(2*np.pi*g.xc[:,None,None]/(g.nx*g.dx))
    base=(s.state.u.copy(),s.state.v.copy(),s.state.w.copy())
    native=momentum_advection_tendency(s.state,g,order=2)
    frozen=frozen_momentum_advection_tendency(base,base,g,order=2)
    eps_rows=[]
    for eps in (1e-2,1e-4,1e-6,1e-8):
        class State: pass
        p=State(); p.u=(1+eps)*base[0]; p.v=(1+eps)*base[1]; p.w=(1+eps)*base[2]
        fd=tuple((a-b)/eps for a,b in zip(momentum_advection_tendency(p,g,order=2),native))
        eps_rows.append({"epsilon":eps,"relative_error_to_2M":sum(rms(a-2*b) for a,b in zip(fd,native))/max(sum(rms(b) for b in native),1e-30),
                         "relative_difference_from_frozen_partition":sum(rms(a-b) for a,b in zip(fd,frozen))/max(sum(rms(a) for a in fd),1e-30)})
    return {"identity":"M((1+eps)U)=(1+eps)^2 M(U); J_M(U)U=2M(U), while exact source partition L_U(U)=M(U)","epsilon_study":eps_rows}


def main():
    out=ROOT/"outputs"/"vorticity_provenance_propagator_validation"; out.mkdir(parents=True,exist_ok=True)
    one=one_dimensional(); eps=epsilon_study(); refine=temporal_refinement(); long=long_cgrid(); jac=full_jacobian_distinction()
    summary={"scope":"instrument validation only; no physical provenance interpretation","one_dimensional":one,
             "directional_finite_difference":eps,"temporal_refinement":refine,"full_jacobian_distinction":jac,"long_memory_reduced_cgrid":long,
             "formulation":"linear frozen-branch partition of the advected velocity operand with the total velocity as common advector; exact additive decomposition, not the full Burgers Jacobian",
             "source_sha256":hashlib.sha256((ROOT/"src/storm_dynamics/vorticity_provenance.py").read_bytes()).hexdigest()}
    (out/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    fig,ax=plt.subplots(figsize=(9,5)); ax.semilogy([r["step"] for r in long["history"]],[r["kappa_omega"] for r in long["history"]],label="$\\kappa_\\omega$")
    ax.set(xlabel="Passo",ylabel="Índice de condicionamento",title="Memória longa reduzida — C-grid"); ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(out/"long_memory_kappa.png",dpi=170); plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,5)); ax.loglog([r["epsilon"] for r in eps],[r["max_abs_error"] for r in eps],"o-"); ax.invert_xaxis(); ax.set(xlabel="$\\epsilon$",ylabel="Erro máximo",title="Derivada direcional do ramo minmod congelado"); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(out/"directional_derivative.png",dpi=170); plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k!="long_memory_reduced_cgrid"},indent=2)); print(json.dumps({k:v for k,v in long.items() if k!="history"},indent=2))

if __name__=="__main__": main()
