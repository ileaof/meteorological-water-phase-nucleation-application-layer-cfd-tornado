"""Tracking must lose identity rather than select an unrelated stronger peak."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('h5py')
spec=importlib.util.spec_from_file_location('diagnostic_analysis',Path(__file__).resolve().parents[1]/'scripts/analyze_diagnostic_sequence.py')
analysis=importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)


def test_tracker_does_not_jump_to_disjoint_peak():
    old=np.zeros((10,10,5),bool); old[2:4,2:4,:]=True
    labels=np.zeros(old.shape,int); labels[7:9,7:9,:]=1
    mask,status,_=analysis.associate(labels,old,None,None,None)
    assert mask is None and status=='lost'
    labels[3:5,2:4,:]=2
    mask,status,_=analysis.associate(labels,old,None,None,None)
    assert status=='overlap' and mask[3,2,0] and not mask[7,7,0]
    split=np.zeros(old.shape,int); split[2,2:4,:]=1; split[3,2:4,:]=2
    mask,status,_=analysis.associate(split,old,None,None,None)
    assert mask is None and status=='ambiguous'


def test_rk4_solid_rotation_on_native_faces(tmp_path):
    import h5py
    from meteorological_flow.grid import Grid
    grid=Grid(12,12,8,1200.,1200.,800.)
    with h5py.File(tmp_path/'analytic.h5','w') as f:
        group=f.create_group('grid')
        for name in ('xc','yc','zc','xf','yf','zf'): group[name]=getattr(grid,name)
        left=f.create_group('a'); right=f.create_group('b')
        omega=.01
        for g,t in ((left,0.),(right,30.)):
            g.attrs['time_s']=t
            g['u']=np.broadcast_to(-omega*(grid.yc[None,:,None]-600),grid.u_shape)
            g['v']=np.broadcast_to(omega*(grid.xc[:,None,None]-600),grid.v_shape)
            g['w']=np.zeros(grid.w_shape)
        interval=analysis.VelocityInterval(left,right,f)
        initial=np.array([[700.,600.,300.]])
        final=analysis.rk4(interval,initial,0.,30.,1.)
        expected=np.array([[600+100*np.cos(.3),600+100*np.sin(.3),300.]])
        assert np.max(abs(final-expected))<1e-7
        returned=analysis.rk4(interval,final,30.,0.,1.)
        assert np.max(abs(returned-initial))<1e-7


def test_end_to_end_synthetic_sequence(tmp_path):
    import h5py
    from meteorological_flow.grid import Grid
    from storm_dynamics.diagnostic_capture import kinematics
    grid=Grid(12,12,12,12000.,12000.,3000.)
    xu,yu,zu=np.meshgrid(grid.xf,grid.yc,grid.zc,indexing='ij')
    xv,yv,zv=np.meshgrid(grid.xc,grid.yf,grid.zc,indexing='ij')
    velocity=(-.008*(yu-6000)*np.exp(-((xu-6000)**2+(yu-6000)**2)/4000**2),
              .008*(xv-6000)*np.exp(-((xv-6000)**2+(yv-6000)**2)/4000**2),
              np.full(grid.w_shape,.1))
    kin=kinematics(velocity,(grid.xc,grid.yc,grid.zc))
    with h5py.File(tmp_path/'synthetic.h5','w') as f:
        f.attrs['budget_nz_with_halo']=grid.nz
        gg=f.create_group('grid')
        for name in ('xc','yc','zc','xf','yf','zf'): gg[name]=getattr(grid,name)
        bg=f.create_group('base'); bg['rho0']=np.ones(grid.nz); bg['rho0_wface']=np.ones(grid.nz+1)
        snapshots=f.create_group('snapshots')
        for i in range(3):
            g=snapshots.create_group(str(i)); g.attrs['time_s']=i*30.
            for name,v in zip(('u','v','w'),velocity): g[name]=v
            for name in ('p','p_dyn','ql','qi'): g[name]=np.zeros(grid.center_shape)
            if i:
                inc=g.create_group('increments')
                for stage in ('advection','surface_drag','les','projection'):
                    group=inc.create_group(stage)
                    for name,v in zip(('u','v','w'),velocity): group[name]=np.zeros_like(v)
                group=g.create_group('kinematic_integrals')
                for name in ('vector_advection','vector_stretching','vector_dilatation','stretching','tilting'): group[name]=kin[name]*30
        summary,peak,seeds,_=analysis.audit_and_track(f,tmp_path,30.)
        assert summary['closure_max_s2']==0
        positions=analysis.trajectories(f,peak,seeds)
        result=analysis.lagrangian_budget(f,positions,tmp_path)
        assert result['valid_budget_intervals']>0
        assert (tmp_path/'pressure_vs_visible_cloud.png').exists()
        summary['metadata']={'config':{'sim':{'grid':{'nx':12,'ny':12,'nz':12}}}}
        summary['lagrangian']={}
        summary['trajectory_sensitivity']={name:0. for name in ('rk4_2s_vs_1s_median_m','rk4_2s_vs_1s_max_m','cadence_30_vs_60s_median_m','cadence_30_vs_60s_max_m')}
        analysis.write_report(summary,result,[],tmp_path,f)
        assert 'comprovado' in (tmp_path/'RELATORIO_FINAL.md').read_text(encoding='utf-8')
