"""Local microphysics only; no causal storm cases are executed."""
from dataclasses import asdict
import numpy as np
import pytest
from precip_microphysics import processes as proc, thermo as th
from precip_microphysics.config import MicrophysicsConfig, from_dict
from precip_microphysics.scheme import BulkMicrophysics
from precip_microphysics.state import MicrophysicsState, SPECIES
from reference_rain_evaporation_legacy import rain_evaporation as legacy


@pytest.fixture(params=['cpu','gpu'])
def xp(request):
    if request.param=='cpu': return np
    cp=pytest.importorskip('cupy')
    try:
        if cp.cuda.runtime.getDeviceCount()<1: pytest.skip('No CUDA GPU')
    except cp.cuda.runtime.CUDARuntimeError: pytest.skip('CUDA unavailable')
    return cp


def state(xp):
    T=xp.asarray([250.,265.,272.,280.,290.,300.])
    P=xp.full(T.shape,85000.)
    return MicrophysicsState(T=T,P=P,rho=xp.ones_like(T),w=xp.full(T.shape,25.),
        qv=.65*th.qsat_water(T,P,xp=xp),qc=xp.full(T.shape,.001),qr=xp.full(T.shape,.002),
        qi=xp.full(T.shape,.001),qs=xp.full(T.shape,.001),qg=xp.full(T.shape,.001),qh=xp.full(T.shape,.001),xp=xp)


def equal_transfers(a,b,xp):
    assert len(a)==len(b)
    for aa,bb in zip(a,b):
        assert (aa.src,aa.dst,aa.name)==(bb.src,bb.dst,bb.name)
        assert bool(xp.array_equal(aa.dq,bb.dq))


@pytest.mark.parametrize('factor',[0.,.5,1.,1.1])
def test_scaling_uncapped_and_legacy_identity(xp,factor):
    st=state(xp); dt=1e-4
    baseline=legacy(st,MicrophysicsConfig(),dt)
    result=proc.rain_evaporation(st,MicrophysicsConfig(rain_evaporation_factor=factor),dt)
    if factor==0:
        assert result==[]
    else:
        xp.testing.assert_allclose(result[0].dq,factor*baseline[0].dq,rtol=3e-16,atol=0)
        if factor==1: equal_transfers(result,baseline,xp)


def test_all_other_rates_invariant_and_active_ice(xp):
    st=state(xp); active_ice=False
    for fn in proc.PROCESS_ORDER:
        if fn is proc.rain_evaporation: continue
        baseline=fn(st,MicrophysicsConfig(),.5)
        if any(tr.src in ('qi','qs','qg','qh') or tr.dst in ('qi','qs','qg','qh') for tr in baseline): active_ice=True
        for factor in (0.,.5,1.1):
            equal_transfers(baseline,fn(st,MicrophysicsConfig(rain_evaporation_factor=factor),.5),xp)
    assert active_ice


def test_full_microphysics_step_bitwise_legacy(xp,monkeypatch):
    modern=state(xp); old=state(xp)
    current=BulkMicrophysics(MicrophysicsConfig()); previous=BulkMicrophysics(MicrophysicsConfig())
    order=proc.PROCESS_ORDER
    oldorder=tuple(legacy if f is proc.rain_evaporation else f for f in order)
    for _ in range(3):
        monkeypatch.setattr(proc,'PROCESS_ORDER',order); current.step(modern,.5)
        monkeypatch.setattr(proc,'PROCESS_ORDER',oldorder); previous.step(old,.5)
        for name in (*SPECIES,'T'):
            assert bool(xp.array_equal(getattr(modern,name),getattr(old,name))),name


def test_caps_remain_enforced(xp):
    st=state(xp)
    for factor in (.5,1.,1.1,100.):
        tr=proc.rain_evaporation(st,MicrophysicsConfig(rain_evaporation_factor=factor),1e6)[0]
        assert bool(xp.all(tr.dq<=st.qr))
        assert bool(xp.all(tr.dq<=xp.maximum(th.qsat_water(st.T,st.P,xp=xp)-st.qv,0)))


@pytest.mark.parametrize('factor',[-1.,float('nan'),float('inf')])
def test_invalid_factors(factor):
    with pytest.raises(ValueError): from_dict({'rain_evaporation_factor':factor})


def test_config_defaults_roundtrip():
    from meteorological_flow.config import from_dict as flow_config
    from meteorological_flow.simulation import _cfg_summary
    assert from_dict({}).rain_evaporation_factor==1.
    assert flow_config({}).physics.rain_evaporation_factor==1.
    assert from_dict(asdict(MicrophysicsConfig(rain_evaporation_factor=.5))).rain_evaporation_factor==.5
    assert flow_config({'physics':{'rain_evaporation_factor':.5}}).physics.rain_evaporation_factor==.5
    assert _cfg_summary(flow_config({'physics':{'rain_evaporation_factor':.5}}))['rain_evaporation_factor']==.5


@pytest.mark.parametrize('lowmem',[False,True])
def test_storm_wires_factor_without_running_storm(lowmem):
    from storm_dynamics.config import build_storm_config
    from storm_dynamics.core import StormSimulation
    n=42 if lowmem else 8; nz=40 if lowmem else 8
    cfg=build_storm_config(nx=n,ny=n,nz=nz,Lx=16000,Ly=16000,Lz=8000,z_stretch=1.05,device='cpu')
    cfg.sim.physics.rain_evaporation_factor=.5
    sim=StormSimulation(cfg)
    assert sim._lowmem_pressure==lowmem
    assert sim.coupler.cfg.rain_evaporation_factor==.5
