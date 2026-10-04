"""Integration tests of the passive observer, distinct from the mature replay."""
import numpy as np
import pytest

from storm_dynamics.config import build_storm_config
from storm_dynamics.core import StormSimulation
from storm_dynamics.muscl_capture import DirectMusclCapture, TERMS, flux_increments
from storm_dynamics.momentum import _u_tendency, _v_tendency


def test_directional_reconstruction_matches_native_tendencies():
    cfg=build_storm_config(nx=8,ny=8,nz=8,Lx=8000.,Ly=8000.,Lz=4000.,z_stretch=1.05,device='cpu')
    sim=StormSimulation(cfg);flux={}
    u=_u_tendency(sim.state,sim.grid,2,sim.grid.periodic,flux)
    v=_v_tendency(sim.state,sim.grid,2,sim.grid.periodic,flux)
    dirs=flux_increments(flux,sim.grid,.125)
    for i,expected in enumerate((u,v)):
        np.testing.assert_allclose(sum(dirs[d][i] for d in 'xyz'),.125*expected,rtol=1e-13,atol=1e-14)


@pytest.mark.parametrize('device',['cpu','gpu'])
def test_observer_is_passive_and_reconstructs_stage_updates(device):
    if device=='gpu':
        cupy=pytest.importorskip('cupy')
        if cupy.cuda.runtime.getDeviceCount()==0:pytest.skip('no GPU')
    cfg=build_storm_config(nx=8,ny=8,nz=8,Lx=8000.,Ly=8000.,Lz=4000.,z_stretch=1.05,device=device)
    plain=StormSimulation(cfg);observed=StormSimulation(cfg)
    ma=np.zeros((8,8));ma[2:5,2:5]=1
    mb=np.zeros((8,8));mb[3:6,2:5]=1
    masks=np.stack((mb,.5*(ma+mb),ma))
    observer=DirectMusclCapture(observed,masks,nk=4)
    observed.diagnostic_observer=observer
    fields=('u','v','w','p','theta','qv','ql','qi','qr','qs','qg','qh')
    for _ in range(3):
        for sim in (plain,observed):
            sim._step(.05);sim.step+=1;sim.t=float(sim.state.t)
        for field in fields:
            assert bool(observed.grid.xp.array_equal(getattr(plain.state,field),getattr(observed.state,field)))
        assert observer.last_record.shape==(len(TERMS),3,4,3)
    stats,balances=observer.finish_block(((ma,mb,False),(ma,mb,True),(ma,ma,False)))
    assert np.isfinite(stats).all() and np.isfinite(balances).all()
    assert all(g['relative_max']<=1e-10 for g in observer.gates.values())
