"""Passive observer invariance and independent discrete/continuum checks."""
import numpy as np
import pytest
from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.diagnostic_capture import DiagnosticCapture, curl, kinematics


def test_analytic_tilting_and_nonuniform_stretching():
    from meteorological_flow.grid import Grid
    g=Grid(8,8,8,800,800,800,z_stretch=1.1)
    xf,yc,zc=np.meshgrid(g.xf,g.yc,g.zc,indexing='ij')
    xc,yf,zc2=np.meshgrid(g.xc,g.yf,g.zc,indexing='ij')
    xc2,yc2,zf=np.meshgrid(g.xc,g.yc,g.zf,indexing='ij')
    native=(-.01*yc-.002*xf+.003*zc,.01*xc-.002*yf+.005*zc2,
            .004*zf+.007*xc2+.011*yc2)
    k=kinematics(native,(g.xc,g.yc,g.zc))
    assert np.allclose(k['omega'][2],.02)
    assert np.allclose(k['stretching'],.00008)
    assert np.allclose(k['tilting'],-.005*.007+.003*.011)


def test_lowmem_observer_preserves_dynamic_and_thermodynamic_pressure(tmp_path):
    pytest.importorskip('h5py')
    cfg=build_storm_config(nx=42,ny=42,nz=40,Lx=16000,Ly=16000,Lz=8000,z_stretch=1.05)
    plain=StormSimulation(cfg); observed=StormSimulation(cfg)
    assert observed._lowmem_pressure
    capture=DiagnosticCapture(observed,tmp_path/'lowmem.h5',{},interval=1.)
    observed.diagnostic_observer=capture
    for sim in (plain,observed):
        sim._step(.1); sim.step+=1; sim.t=float(sim.state.t)
    for name in ('u','v','w','p','p_dyn','theta','qv','ql','P_total'):
        assert np.array_equal(getattr(plain.state,name),getattr(observed.state,name)),name
    capture.close(observed)


def test_continuous_material_terms_are_not_flux_form_in_divergent_flow():
    """A nonzero remainder need not mean dissipation: flux adds -curl(u div u)."""
    from meteorological_flow.grid import Grid
    from meteorological_flow.state import FlowState
    from storm_dynamics.momentum import momentum_advection_tendency
    g=Grid(20,20,20,2000,2000,2000)
    st=FlowState.zeros(g)
    st.u[:]=-.01*(g.yc[None,:,None]-1000)
    st.v[:]=.01*(g.xc[:,None,None]-1000)
    st.w[:]=.001*g.zf[None,None,:]
    velocity=(st.u,st.v,st.w); coords=(g.xc,g.yc,g.zc)
    kin=kinematics(velocity,coords)
    actual=curl(momentum_advection_tendency(st,g,order=2,periodic=False),coords)[2]
    material=(kin['vector_advection']+kin['vector_stretching']+kin['vector_dilatation'])[2]
    region=np.s_[4:-4,4:-4,4:-4]
    assert np.allclose(material[region],0,atol=1e-15)
    assert np.allclose(actual[region],-.02*.001,atol=1e-15)


@pytest.mark.parametrize('device',['cpu','gpu'])
def test_observer_is_bitwise_passive_and_budget_telescopes(tmp_path,device):
    h5py=pytest.importorskip('h5py')
    if device=='gpu':
        cp=pytest.importorskip('cupy')
        try: cp.zeros(1)
        except Exception as exc: pytest.skip(str(exc))
    cfg=build_storm_config(nx=10,ny=10,nz=8,Lx=8000,Ly=8000,Lz=4000,z_stretch=1.05,device=device)
    plain=StormSimulation(cfg); observed=StormSimulation(cfg)
    capture=DiagnosticCapture(observed,tmp_path/'sequence.h5',{},interval=.3)
    observed.diagnostic_observer=capture
    for i in range(4):
        for sim in (plain,observed):
            sim._step(.1); sim.step+=1; sim.t=float(sim.state.t)
        for name in ('u','v','w','p','theta','qv','ql','qi','qr'):
            assert np.array_equal(plain.backend.to_cpu(getattr(plain.state,name)),observed.backend.to_cpu(getattr(observed.state,name))),name
    capture.close(observed)
    with h5py.File(tmp_path/'sequence.h5','r') as f:
        nk=int(f.attrs['budget_nz_with_halo']); coord=tuple(f['grid'][n][:] for n in ('xc','yc','zc'))
        coord=(*coord[:2],coord[2][:nk]); snapshots=list(f['snapshots'].values())
        for a,b in zip(snapshots,snapshots[1:]):
            vel=lambda group: tuple(group[n][:,:,:nk+(n=='w')] for n in ('u','v','w'))
            tendency=sum((curl(vel(stage),coord) for stage in b['increments'].values()))
            delta=curl(vel(b),coord)-curl(vel(a),coord)
            assert np.max(abs(delta-tendency))<1e-14
        assert f.attrs['status']=='complete'
        assert len(f['steps'])==4
        # .3 is saved inside the stage callback, .4 only by close().
        assert [int(g.attrs['step']) for g in snapshots] == [0, 3, 4]
        assert snapshots[-1].attrs['time_s'] == observed.state.t
        assert snapshots[-1].attrs['step'] == int(f['steps'][-1,2])+1
        assert 'increments' in snapshots[-1]
