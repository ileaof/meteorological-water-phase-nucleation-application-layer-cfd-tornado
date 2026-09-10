"""Read-only saved-field audit; no solver imports or simulation integration.

Run from repository root. Outputs go to outputs/implementation_audit/.
Derivatives use actual coordinates, second-order one-sided edges; statistics
exclude lateral outer 1/6 (a sensitivity mask, not a proven sponge boundary).
"""
from pathlib import Path
import csv
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/implementation_audit'


def derivative(a, coordinate, axis):
    return np.gradient(a, coordinate, axis=axis, edge_order=2)


def fields(u, v, w, x, y, z):
    ux, uy, uz = [derivative(u, c, i) for i, c in enumerate((x,y,z))]
    vx, vy, vz = [derivative(v, c, i) for i, c in enumerate((x,y,z))]
    wx, wy, wz = [derivative(w, c, i) for i, c in enumerate((x,y,z))]
    xi, eta, zeta = wy-vz, uz-wx, vx-uy
    convergence = -(ux+vy)
    return dict(xi=xi, eta=eta, zeta=zeta, convergence=convergence,
                stretching=zeta*wz, tilting=xi*wx+eta*wy,
                dilatation=-zeta*(ux+vy+wz), net_concentration=zeta*convergence)


def self_test():
    x=np.arange(7.)*100; y=np.arange(8.)*120; z=np.array([10.,30.,80.,160.,300.])
    X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
    a=fields(-.01*Y-.002*X, .01*X-.002*Y, .004*Z, x,y,z)
    assert np.allclose(a['zeta'],.02)
    assert np.allclose(a['convergence'],.004)
    assert np.allclose(a['stretching'],.00008)
    assert np.allclose(a['tilting'],0)
    b=fields(.003*Z, .005*Z, .007*X+.011*Y,x,y,z)
    assert np.allclose(b['tilting'], -.005*.007+.003*.011)


def audit(path):
    with np.load(path,allow_pickle=False) as data:
        if not all(k in data for k in ('u','v','w')): return None
        f={k:data[k] for k in data.files}
    u,v,w=[np.asarray(f[k],float) for k in ('u','v','w')]
    shape=f['theta'].shape if 'theta' in f else w.shape
    if u.shape[0]==shape[0]+1: u=(u[:-1]+u[1:])/2
    if v.shape[1]==shape[1]+1: v=(v[:,:-1]+v[:,1:])/2
    if w.shape[2]==shape[2]+1: w=(w[:,:,:-1]+w[:,:,1:])/2
    assert u.shape==v.shape==w.shape, (path,u.shape,v.shape,w.shape)
    nx,ny,nz=u.shape
    if all(k in f for k in ('xc','yc','zc')):
        x,y,z=[np.asarray(f[k],float) for k in ('xc','yc','zc')]
        ratio=float(f.get('dx',np.diff(x).mean()))/np.diff(x).mean()
        scale=1000 if np.isclose(ratio,1000) else 1
        x,y,z=[a*scale for a in (x,y,z)]
        provenance=f'saved coordinates; scale to metres={scale}'
    elif 'parent' in path.name:
        # Confirmed in cache-producing scripts and render_storm_3d.py.
        x=(np.arange(nx)+.5)*600; y=(np.arange(ny)+.5)*600
        edges=np.r_[0,np.cumsum(1.05**np.arange(nz))]; edges*=15000/edges[-1]
        z=(edges[:-1]+edges[1:])/2
        provenance='reconstructed parent: dx=dy=600m, Lz=15000m, stretch=1.05; cache lacks grid metadata'
    else: raise ValueError(f'No verified grid: {path}')
    a=fields(u,v,w,x,y,z)
    bx,by=max(2,nx//6),max(2,ny//6)
    interior=np.zeros((nx,ny),bool); interior[bx:-bx,by:-by]=True
    rows=[]; previous=None
    for k,zk in enumerate(z):
        zz=a['zeta'][:,:,k]; ww=w[:,:,k]
        ij=np.unravel_index(np.argmax(np.where(interior,zz,-np.inf)),zz.shape)
        iw=np.unravel_index(np.argmax(np.where(interior,ww,-np.inf)),ww.shape)
        positive=np.maximum(zz[interior],0); up=np.maximum(ww[interior],0)
        denom=np.linalg.norm(positive)*np.linalg.norm(up)
        row=dict(z_m=float(zk),x_zeta_m=float(x[ij[0]]),y_zeta_m=float(y[ij[1]]),
                 zeta_max_s=float(zz[ij]),w_at_zeta_ms=float(ww[ij]),w_max_ms=float(ww[iw]),
                 peak_offset_m=float(np.hypot(x[ij[0]]-x[iw[0]],y[ij[1]]-y[iw[1]])),
                 positive_w_zeta_cosine=float(np.dot(positive,up)/denom) if denom else 0.,
                 adjacent_peak_displacement_m=None if previous is None else float(np.hypot(x[ij[0]]-previous[0],y[ij[1]]-previous[1])))
        previous=(x[ij[0]],y[ij[1]])
        for name,arr in a.items():
            row[name+'_at_zeta']=float(arr[:,:,k][ij])
            row[name+'_absmean']=float(np.abs(arr[:,:,k][interior]).mean())
        for name in ('p','p_dyn'):
            if name in f and f[name].shape==u.shape:
                pp=f[name][:,:,k]
                row[name+'_core_minus_plane_median']=float(pp[ij]-np.median(pp[interior]))
        rows.append(row)
    # Follow a local cyclonic ridge down from 1.5 km, with a fixed 1 km
    # horizontal search radius. This is a geometric candidate, not a vortex ID.
    X,Y=np.meshgrid(x,y,indexing='ij')
    anchor=int(np.argmin(abs(z-1500)))
    center=(rows[anchor]['x_zeta_m'],rows[anchor]['y_zeta_m'])
    ridge=[]
    for k in range(anchor,-1,-1):
        mask=interior & (np.hypot(X-center[0],Y-center[1])<=1000)
        ij=np.unravel_index(np.argmax(np.where(mask,a['zeta'][:,:,k],-np.inf)),interior.shape)
        center=(float(x[ij[0]]),float(y[ij[1]]))
        ridge.append(dict(z_m=float(z[k]),x_m=center[0],y_m=center[1],
                          zeta_s=float(a['zeta'][ij[0],ij[1],k]),w_ms=float(w[ij[0],ij[1],k])))
    low=z<=2000
    report=dict(path=str(path.relative_to(ROOT)),keys=list(f),shape=list(u.shape),
                coordinates=provenance,dx_m=float(np.diff(x).mean()),lowest_z_m=float(z[0]),
                levels_0_2km=int(low.sum()),time_s=float(np.ravel(f['_t'])[0]) if '_t' in f else None,
                low_abs_zeta_all=float(np.abs(a['zeta'][:,:,low]).max()),
                low_abs_zeta_interior=float(np.abs(a['zeta'][interior][:,low]).max()),
                rows=rows,downward_ridge_1km_search=ridge)
    report['border_sensitivity_low_abs_zeta']={}
    for fraction in (.05,.10,.20,.25):
        mx,my=max(2,int(nx*fraction)),max(2,int(ny*fraction))
        report['border_sensitivity_low_abs_zeta'][str(fraction)]=float(np.abs(a['zeta'][mx:-mx,my:-my,low]).max())
    for name in ('p','p_dyn'):
        if name in f: report[name+'_range']=[float(f[name].min()),float(f[name].max())]
    if 'zeta' in f and f['zeta'].shape==u.shape:
        report['stored_zeta_interior_max_error']=float(np.abs(f['zeta'][interior]-a['zeta'][interior]).max())
    if 'ql' in f:
        cond=f['ql']+f.get('qi',0)
        report['lowest_condensate_level_m']=float(z[np.where((cond>1e-5).any(axis=(0,1)))[0][0]]) if (cond>1e-5).any() else None
        report['lowest_level_cloud_fraction']=float((cond[:,:,0]>1e-5).mean())
    stem=str(path.relative_to(ROOT/'outputs')).replace('/','__').replace('\\','__').replace('.npz','')
    with (OUT/(stem+'.csv')).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    if path.name=='parent_matured_120_48_2800.npz':
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,3,figsize=(12,5),constrained_layout=True)
        axs[0].plot([r['zeta_max_s'] for r in rows],z/1000); axs[0].set_xlabel('Interior maximum zeta (1/s)')
        axs[1].plot([r['w_at_zeta_ms'] for r in rows],z/1000,label='w at zeta peak')
        axs[1].plot([r['w_max_ms'] for r in rows],z/1000,label='maximum w'); axs[1].legend(); axs[1].set_xlabel('w (m/s)')
        axs[2].plot([r['peak_offset_m']/1000 for r in rows],z/1000); axs[2].set_xlabel('Updraft / zeta peak offset (km)')
        for ax in axs: ax.set_ylim(0,2); ax.set_ylabel('Height (km)'); ax.grid(alpha=.3)
        fig.suptitle('Saved mature parent: independent 0–2 km audit')
        fig.savefig(OUT/'parent_vertical_profiles.png',dpi=160); plt.close(fig)
    return report


if __name__=='__main__':
    self_test(); OUT.mkdir(parents=True,exist_ok=True)
    reports=[]; errors=[]
    for path in sorted((ROOT/'outputs').rglob('*.npz')):
        try:
            result=audit(path)
            if result:
                reports.append(result)
                print(result['path'],result['dx_m'],result['low_abs_zeta_interior'],flush=True)
        except Exception as e: errors.append(dict(path=str(path),error=str(e)))
    (OUT/'audit.json').write_text(json.dumps(dict(analytic_checks='passed',reports=reports,errors=errors),indent=2))
    history=ROOT/'outputs/supercell_alignment/summary.json'
    if history.exists():
        records=json.loads(history.read_text())['history']
        with (OUT/'legacy_alignment_history.csv').open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    print('Audited',len(reports),'errors',errors)
