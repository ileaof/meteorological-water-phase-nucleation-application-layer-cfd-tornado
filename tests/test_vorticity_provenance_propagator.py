"""Extended stability/consistency tests for the frozen-MUSCL provenance map."""
from types import SimpleNamespace

import numpy as np
import pytest

from meteorological_flow.advection import _minmod
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.momentum import momentum_advection_tendency
from storm_dynamics.vorticity_provenance import (
    SOURCE_NAMES,
    VorticityProvenanceTracer,
    _frozen_slope,
    curl_native,
    frozen_momentum_advection_tendency,
    frozen_muscl_advection_tendency,
)


def small_sim(nx=16, ny=12, nz=8):
    cfg = build_storm_config(nx=nx, ny=ny, nz=nz, Lx=16000., Ly=12000.,
                             Lz=5000., z_stretch=1., duration=1., device="cpu")
    return StormSimulation(cfg)


@pytest.mark.parametrize("kind", ["smooth", "pulse", "gradient", "sign", "limited"])
def test_1d_many_crossings_is_stable_and_additive(kind):
    sim=small_sim(nx=64,ny=4,nz=4); g=sim.grid
    x=np.arange(g.nx,dtype=float)
    if kind=="smooth": q=np.sin(2*np.pi*x/g.nx)
    elif kind=="pulse": q=((x>=18)&(x<34)).astype(float)
    elif kind=="gradient": q=(x-g.nx/2)/g.nx
    elif kind=="sign": q=np.where(x<g.nx/2,-1.,1.)
    else: q=np.sin(2*np.pi*x/g.nx)+.35*np.sin(10*np.pi*x/g.nx)
    q=np.broadcast_to(q[:,None,None],(g.nx,g.ny,g.nz)).copy()
    total=q.copy(); a=.55*q+.08*np.roll(q,3,axis=0); b=total-a
    vel=(np.full(g.u_shape,12.),np.zeros(g.v_shape),np.zeros(g.w_shape))
    dt=.35*g.dx/12.; initial=max(np.max(abs(total)),1.)
    # >8 domain crossing times.
    for _ in range(1600):
        tn=frozen_muscl_advection_tendency(total,total,vel,g)
        an=frozen_muscl_advection_tendency(a,total,vel,g)
        bn=frozen_muscl_advection_tendency(b,total,vel,g)
        total+=dt*tn; a+=dt*an; b+=dt*bn
    assert np.max(np.abs(total-(a+b))) < 2e-12
    assert np.all(np.isfinite(total))
    assert np.max(np.abs(total)) <= 1.05*initial


def test_old_centered_euler_reference_amplifies_fourier_mode():
    n=64; c=.45; x=np.arange(n); q=np.sin(2*np.pi*9*x/n); q0=np.linalg.norm(q)
    for _ in range(300):
        q -= .5*c*(np.roll(q,-1)-np.roll(q,1))
    assert np.linalg.norm(q) > 100*q0


@pytest.mark.parametrize("flow", ["translation", "shear", "rotation", "deformation"])
def test_multidimensional_passive_transport_is_additive_without_explosion(flow):
    sim=small_sim(nx=20,ny=18,nz=10); g=sim.grid
    x=g.xc[:,None,None]; y=g.yc[None,:,None]; z=g.zc[None,None,:]
    q=np.sin(2*np.pi*x/(g.nx*g.dx))*np.cos(2*np.pi*y/(g.ny*g.dy))*np.exp(-z/6000.)
    u=np.zeros(g.u_shape); v=np.zeros(g.v_shape); w=np.zeros(g.w_shape)
    if flow=="translation": u[:]=7.; v[:]=3.; w[:,:,1:-1]=.4
    elif flow=="shear": u[:]=4.+5*(g.yc[None,:,None]/(g.ny*g.dy)-.5)
    elif flow=="rotation":
        u[:]=-(g.yc[None,:,None]-g.ny*g.dy/2)/3000.
        v[:]=(g.xc[:,None,None]-g.nx*g.dx/2)/3000.
    else:
        # Periodic, analytically divergence-free strain field.
        lx=g.nx*g.dx; ly=g.ny*g.dy
        u[:]=3*np.sin(2*np.pi*g.xf[:,None,None]/lx)*np.cos(2*np.pi*g.yc[None,:,None]/ly)
        v[:]=-3*(ly/lx)*np.cos(2*np.pi*g.xc[:,None,None]/lx)*np.sin(2*np.pi*g.yf[None,:,None]/ly)
    speed=max(np.max(abs(u)),np.max(abs(v)),np.max(abs(w)),1e-9)
    dt=.12*min(g.dx,g.dy,g.dz)/speed
    total=q.copy(); a=.4*q+.03*np.roll(q,1,0); b=total-a; initial=np.linalg.norm(total)
    for _ in range(250):
        tn=frozen_muscl_advection_tendency(total,total,(u,v,w),g)
        an=frozen_muscl_advection_tendency(a,total,(u,v,w),g)
        bn=frozen_muscl_advection_tendency(b,total,(u,v,w),g)
        total+=dt*tn; a+=dt*an; b+=dt*bn
    assert np.linalg.norm(total-(a+b)) <= 2e-12*max(initial,1.)
    assert np.linalg.norm(total) < 3*initial


def test_cgrid_curl_and_native_map_close_for_random_split():
    sim=small_sim(); g=sim.grid; total=(sim.state.u.copy(),sim.state.v.copy(),sim.state.w.copy())
    rng=np.random.default_rng(20260908)
    parts=[]
    for weight in (.2,.3): parts.append(tuple(weight*q+.01*rng.normal(size=q.shape) for q in total))
    parts.append(tuple(q-parts[0][i]-parts[1][i] for i,q in enumerate(total)))
    native=momentum_advection_tendency(sim.state,g,order=2)
    split=[frozen_momentum_advection_tendency(p,total,g,order=2) for p in parts]
    summed=tuple(sum(s[i] for s in split) for i in range(3))
    for n,s in zip(native,summed): assert np.allclose(n,s,rtol=3e-13,atol=3e-13)
    curl_sum=sum(curl_native(p,g) for p in parts)
    assert np.allclose(curl_native(total,g),curl_sum,rtol=2e-13,atol=2e-13)


def test_frozen_minmod_matches_directional_derivative_away_from_switches():
    sim=small_sim(nx=32,ny=4,nz=4); g=sim.grid
    x=np.arange(g.nx,dtype=float)[:,None,None]
    total=(x+.02*x*x)*np.ones((1,g.ny,g.nz)); direction=np.sin(.17*x)*np.ones((1,g.ny,g.nz))
    analytic=_frozen_slope(direction,total,0,g,True)
    errors=[]
    base=_minmod(total-np.roll(total,1,0),np.roll(total,-1,0)-total)
    for eps in (1e-2,1e-4,1e-6,1e-8):
        pert=total+eps*direction
        fd=(_minmod(pert-np.roll(pert,1,0),np.roll(pert,-1,0)-pert)-base)/eps
        errors.append(np.max(abs(fd-analytic)))
    assert min(errors[:3]) < 2e-8
    assert errors[-1] < 2e-6


def test_full_burgers_jacobian_is_distinct_from_exact_frozen_advector_partition():
    sim=small_sim(); g=sim.grid
    sim.state.u += .8*np.sin(2*np.pi*g.xf[:,None,None]/(g.nx*g.dx))*np.cos(2*np.pi*g.yc[None,:,None]/(g.ny*g.dy))
    sim.state.v += -.8*(g.ny*g.dy/(g.nx*g.dx))*np.cos(2*np.pi*g.xc[:,None,None]/(g.nx*g.dx))*np.sin(2*np.pi*g.yf[None,:,None]/(g.ny*g.dy))
    sim.state.w[:,:,1:-1] += .1*np.sin(2*np.pi*g.xc[:,None,None]/(g.nx*g.dx))
    base=(sim.state.u.copy(),sim.state.v.copy(),sim.state.w.copy())
    native=momentum_advection_tendency(sim.state,g,order=2)
    frozen=frozen_momentum_advection_tendency(base,base,g,order=2)
    for a,b in zip(native,frozen): assert np.array_equal(a,b)
    # M((1+eps)U)=(1+eps)^2 M(U) while L_U(U)=M(U): the full
    # directional derivative is 2M(U), so it cannot be used as an additive
    # source partition whose labels sum back to M(U).
    eps=1e-6; pert=SimpleNamespace(u=(1+eps)*base[0],v=(1+eps)*base[1],w=(1+eps)*base[2])
    fd=tuple((p-n)/eps for p,n in zip(momentum_advection_tendency(pert,g,order=2),native))
    relative=sum(np.sqrt(np.mean((d-2*n)**2)) for d,n in zip(fd,native))/sum(np.sqrt(np.mean(n*n)) for n in native)
    assert relative < 2e-6
    assert sum(np.sqrt(np.mean((d-f)**2)) for d,f in zip(fd,frozen)) > .45*sum(np.sqrt(np.mean(d*d)) for d in fd)


@pytest.mark.parametrize("stage,source", [
    ("buoyancy","buoyancy"),("les","les"),("surface_drag","surface_drag"),
    ("coriolis","coriolis"),("projection","projection")])
def test_isolated_native_increment_is_born_only_in_its_label_and_is_transportable(stage,source):
    sim=small_sim(); tracer=VorticityProvenanceTracer(sim); rng=np.random.default_rng(17)
    delta=tuple(rng.normal(scale=1e-5,size=q.shape) for q in (sim.state.u,sim.state.v,sim.state.w))
    delta[2][:,:,0]=0.; delta[2][:,:,-1]=0.
    for name,d in zip(("u","v","w"),delta): setattr(sim.state,name,getattr(sim.state,name)+d)
    tracer.mark(sim,stage,.1)
    for name in SOURCE_NAMES:
        if name in ("initial",source): continue
        assert all(np.count_nonzero(q[tracer.source_index[name]])==0 for q in tracer.velocity_sources)
    before=tuple(q[tracer.source_index[source]].copy() for q in tracer.velocity_sources)
    from storm_dynamics.momentum import add_momentum_advection
    add_momentum_advection(sim.state,sim.grid,.02,order=2)
    tracer.mark(sim,"advection",.02)
    after=tuple(q[tracer.source_index[source]] for q in tracer.velocity_sources)
    assert any(np.max(abs(a-b))>0 for a,b in zip(after,before))
    m=tracer._closure_metrics(); assert m["omega_relative_rms"]<1e-11


def test_long_memory_reduced_cgrid_has_no_secular_conditioning_growth():
    sim=small_sim(nx=10,ny=10,nz=6); g=sim.grid
    rng=np.random.default_rng(9)
    total=(sim.state.u.copy(),sim.state.v.copy(),sim.state.w.copy())
    parts=[tuple(.5*q+.001*rng.normal(size=q.shape) for q in total)]
    parts.append(tuple(q-parts[0][i] for i,q in enumerate(total)))
    kappas=[]; dt=.005
    for n in range(2500):
        tend_total=frozen_momentum_advection_tendency(total,total,g,order=2)
        tend_parts=[frozen_momentum_advection_tendency(p,total,g,order=2) for p in parts]
        total=tuple(q+dt*d for q,d in zip(total,tend_total))
        parts=[tuple(q+dt*d for q,d in zip(p,t)) for p,t in zip(parts,tend_parts)]
        if n%100==0:
            o=curl_native(total,g); os=[curl_native(p,g) for p in parts]
            scale=np.sqrt(np.mean(o*o)); kappas.append(sum(np.sqrt(np.mean(x*x)) for x in os)/max(scale,1e-30))
            assert np.max(abs(o-sum(os))) < 2e-12
    assert np.all(np.isfinite(kappas))
    assert max(kappas) < 3*max(kappas[0],1.)
