"""Offline geometry/source audit of the existing 27 parcels; no solver calls.

Reads the immutable fixed-grid sequence and previous parcel IDs/positions.
No trajectories are reseeded and no simulation is integrated.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import numpy as np
import h5py
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from storm_dynamics.diagnostic_capture import centered, grad, curl


def read_csv(path):
    with path.open(encoding='utf-8-sig') as stream: return list(csv.DictReader(stream))


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        writer.writeheader(); writer.writerows(rows)


def wrap(angle): return (angle+np.pi)%(2*np.pi)-np.pi


def angular_decomposition(omega0,omega1,g0,g1):
    """Exact symmetric decomposition of delta(omega.g), including magnitudes.

    Orientation terms average the two orders of freezing/rotating vectors
    (Shapley split of cos(beta-alpha)); these are kinematic counterfactuals,
    not model interventions or causal attribution to a force.
    """
    o0=np.linalg.norm(omega0); o1=np.linalg.norm(omega1)
    h0=np.linalg.norm(g0); h1=np.linalg.norm(g1)
    if min(o0,o1,h0,h1)<1e-10: return None
    a0=np.arctan2(omega0[1],omega0[0]); a1=np.arctan2(omega1[1],omega1[0])
    b0=np.arctan2(g0[1],g0[0]); b1=np.arctan2(g1[1],g1[0])
    c00=np.cos(b0-a0); c10=np.cos(b0-a1); c01=np.cos(b1-a0); c11=np.cos(b1-a1)
    orientation_omega=.5*((c10-c00)+(c11-c01))
    orientation_gradient=.5*((c01-c00)+(c11-c10))
    pmid=.5*(o0*h0+o1*h1); cmid=.5*(c00+c11)
    parts=dict(omega_orientation=pmid*orientation_omega,gradient_orientation=pmid*orientation_gradient,
               omega_magnitude=cmid*(o1-o0)*(h0+h1)/2,
               gradient_magnitude=cmid*(h1-h0)*(o0+o1)/2)
    delta=float(np.dot(omega1,g1)-np.dot(omega0,g0))
    return dict(**parts,delta_T=delta,closure=delta-sum(parts.values()),
                delta_alpha_deg=float(np.degrees(wrap(a1-a0))),delta_beta_deg=float(np.degrees(wrap(b1-b0))),
                cos0=c00,cos1=c11,
                cos_if_only_omega_rotated=c10,cos_if_only_gradient_rotated=c01)


def trilinear_gradient(a,coords,points):
    """Analytic derivative of the cell-centred trilinear interpolant.

    Alternate estimator for sensitivity, not a higher-resolution solution.
    """
    indices=[]; weights=[]; spans=[]
    for axis,c in enumerate(coords):
        idx=np.clip(np.searchsorted(c,points[:,axis],side='right')-1,0,len(c)-2)
        indices.append(idx); spans.append(c[idx+1]-c[idx]); weights.append((points[:,axis]-c[idx])/spans[-1])
    result=np.zeros((len(points),3))
    for ix in (0,1):
        for iy in (0,1):
            for iz in (0,1):
                corner=a[indices[0]+ix,indices[1]+iy,indices[2]+iz]
                bits=(ix,iy,iz)
                for axis in range(3):
                    factor=np.ones(len(points))*(1 if bits[axis] else -1)/spans[axis]
                    for other in range(3):
                        if other!=axis: factor*=weights[other] if bits[other] else 1-weights[other]
                    result[:,axis]+=corner*factor
    return result


def field_bundle(group,f,nk):
    coords=tuple(f['grid'][name][:] for name in ('xc','yc','zc')); coords=(*coords[:2],coords[2][:nk])
    velocity=tuple(group[n][:,:,:nk+(n=='w')] for n in ('u','v','w'))
    u,v,w=centered(velocity); ux,uy,uz=grad(u,coords); vx,vy,vz=grad(v,coords); wx,wy,wz=grad(w,coords)
    xi,eta,zeta=wy-vz,uz-wx,vx-uy
    theta=group['theta'][:,:,:nk]; qv=group['qv'][:,:,:nk]
    loading=sum(group[n][:,:,:nk] for n in ('ql','qi','qr','qs','qg','qh'))
    t0=f['base/theta0'][:nk]; q0=f['base/qv0'][:nk]
    metadata=json.loads(f.attrs['metadata']); gravity=metadata['config']['sim']['flow']['gravity']
    moist=metadata['config']['sim']['physics']['moisture_buoyancy']
    bthermal=gravity*(theta-t0)/t0
    bmoist=gravity*.61*(qv-q0) if moist else np.zeros_like(theta)
    bloading=-gravity*loading if moist else np.zeros_like(theta)
    B=bthermal+bmoist+bloading; bx,by,bz=grad(B,coords)
    tv=theta*(1+.61*qv-loading); tvbase=t0*(1+.61*q0)
    Bface=np.empty((*B.shape[:2],nk+1)); Bface[:,:,1:-1]=(B[:,:,:-1]+B[:,:,1:])/2
    Bface[:,:,0]=B[:,:,0]; Bface[:,:,-1]=B[:,:,-1]
    bc=(Bface[:,:,:-1]+Bface[:,:,1:])/2; bcx,bcy,_=grad(bc,coords)
    out=dict(u=u,v=v,w=w,xi=xi,eta=eta,zeta=zeta,wx=wx,wy=wy,
             T=xi*wx+eta*wy,Tx=xi*wx,Ty=eta*wy,convergence=-(ux+vy),
             B=B,Bthermal=bthermal,Bmoist=bmoist,Bloading=bloading,
             Bx=bx,By=by,buoyancy_curl_xi=bcy,buoyancy_curl_eta=-bcx,
             dtheta=theta-t0,dtv_base=tv-tvbase,dtv_mean=tv-tv.mean(axis=(0,1),keepdims=True))
    for name,g in (('wx',wx),('wy',wy)):
        derivatives=grad(g,coords)
        out['transport_'+name]=sum(a*b for a,b in zip((u,v,w),derivatives))
        out['transport_horizontal_'+name]=u*derivatives[0]+v*derivatives[1]
        out['transport_vertical_'+name]=w*derivatives[2]
    return out,coords


def sample(fields,coords,points):
    return {name:RegularGridInterpolator(coords,a,bounds_error=False,fill_value=np.nan)(points) for name,a in fields.items()}


def collect(sequence,prior,out):
    old=read_csv(prior/'parcels.csv'); old_budget=read_csv(prior/'lagrangian_budget.csv')
    points_by_index={}
    for row in old: points_by_index.setdefault(int(row['index']),[]).append(row)
    for rows in points_by_index.values(): rows.sort(key=lambda r:int(r['parcel']))
    existing_budget={(int(r['index']),int(r['parcel'])):r for r in old_budget}
    tracks={int(r['index']):r for r in read_csv(prior/'vortex_track.csv')}
    rows=[]; rates=[]; surfaces=[]; saved_maps={}; cache={}
    with h5py.File(sequence,'r') as f:
        if f.attrs['status']!='complete': raise RuntimeError('Incomplete input sequence')
        nk=int(f.attrs['budget_nz_with_halo']); groups=list(f['snapshots'].values())
        selected=[i for i,g in enumerate(groups) if 2370<=g.attrs['time_s']<=3331]
        map_ids={min(selected,key=lambda i:abs(groups[i].attrs['time_s']-t)) for t in (2490,2790,2940,3090,3300)}
        previous=None
        for i in selected:
            group=groups[i]; time=float(group.attrs['time_s'])
            points=np.array([[float(r[k]) for k in ('x_m','y_m','z_m')] for r in points_by_index[i]])
            fields,coords=field_bundle(group,f,nk); sampled=sample(fields,coords,points)
            ksurface=int(np.argmin(abs(coords[2]-100)))
            cold=fields['dtv_base'][:,:,ksurface]<-1
            dx=float(np.diff(coords[0]).mean()); dy=float(np.diff(coords[1]).mean())
            distance=(distance_transform_edt(~cold,sampling=(dx,dy))-distance_transform_edt(cold,sampling=(dx,dy))) if cold.any() else np.full(cold.shape,np.nan)
            interp2=lambda array:RegularGridInterpolator(coords[:2],array,bounds_error=False,fill_value=np.nan)(points[:,:2])
            distance_p=interp2(distance); groundtv=interp2(fields['dtv_base'][:,:,ksurface])
            groundw=interp2(fields['w'][:,:,ksurface]); groundconv=interp2(fields['convergence'][:,:,ksurface])
            altu=trilinear_gradient(fields['u'],coords,points); altv=trilinear_gradient(fields['v'],coords,points); altw=trilinear_gradient(fields['w'],coords,points)
            altT=(altw[:,1]-altv[:,2])*altw[:,0]+(altu[:,2]-altw[:,0])*altw[:,1]
            surfaces.append(dict(time_s=time,z_surface_sample_m=float(coords[2][ksurface]),cold_area_km2=float(cold.sum()*dx*dy/1e6),
                                 cold_min_dtv_K=float(fields['dtv_base'][:,:,ksurface].min()),
                                 cold_mean_reference_area_km2=float((fields['dtv_mean'][:,:,ksurface]<-1).sum()*dx*dy/1e6),
                                 max_buoyancy_gradient_s2=float(np.hypot(fields['Bx'][:,:,ksurface],fields['By'][:,:,ksurface]).max())))
            for n,point in enumerate(points):
                values={key:float(a[n]) for key,a in sampled.items()}
                omega=np.array([values['xi'],values['eta']]); gg=np.array([values['wx'],values['wy']])
                norms=np.linalg.norm(omega)*np.linalg.norm(gg)
                cosine=float(np.dot(omega,gg)/norms) if norms>1e-15 else np.nan
                rows.append(dict(index=i,parcel=n,time_s=time,x_m=point[0],y_m=point[1],z_m=point[2],**values,
                                 T_product=float(np.dot(omega,gg)),Tx_product=float(omega[0]*gg[0]),Ty_product=float(omega[1]*gg[1]),
                                 T_trilinear_derivative=float(altT[n]),angle_deg=float(np.degrees(np.arccos(np.clip(cosine,-1,1)))),
                                 alpha_deg=float(np.degrees(np.arctan2(omega[1],omega[0]))),beta_deg=float(np.degrees(np.arctan2(gg[1],gg[0]))),
                                 omega_h_s=float(np.linalg.norm(omega)),grad_w_s=float(np.linalg.norm(gg)),
                                 cold_boundary_distance_m=float(distance_p[n]),surface_dtv_K=float(groundtv[n]),
                                 surface_w_ms=float(groundw[n]),surface_convergence_s=float(groundconv[n])))
            cache[i]=dict(time=time,points=points,sampled=sampled)
            if previous is not None:
                pi,pfields,ppoints=previous; dt=time-float(groups[pi].attrs['time_s']); mid=(points+ppoints)/2
                contributions={}
                for name,stage in group['increments'].items():
                    dwf=stage['w'][:]; dwc=(dwf[:,:,:-1]+dwf[:,:,1:])/2
                    gx,gy,_=grad(dwc,coords)
                    contributions[name]=np.column_stack(list(sample({'x':gx,'y':gy},coords,mid).values()))
                for direction in ('horizontal','vertical'):
                    names=('transport_'+direction+'_wx','transport_'+direction+'_wy')
                    contributions['parcel_transport_'+direction]=np.column_stack([.5*dt*(sample({name:fields[name]},coords,mid)[name]+sample({name:pfields[name]},coords,mid)[name]) for name in names])
                g0=np.column_stack([cache[pi]['sampled'][key] for key in ('wx','wy')]); g1=np.column_stack([sampled[key] for key in ('wx','wy')])
                contributions['sampling_residual']=g1-g0-sum(contributions.values())
                for n in range(len(points)):
                    omega0=np.array([cache[pi]['sampled'][key][n] for key in ('xi','eta')]); omega1=np.array([sampled[key][n] for key in ('xi','eta')])
                    omid=(omega0+omega1)/2; gmid=(g0[n]+g1[n])/2
                    row=dict(index=i,parcel=n,time_s=time,dt_s=dt,
                             omega_angle_change_rad=float(wrap(np.arctan2(omega1[1],omega1[0])-np.arctan2(omega0[1],omega0[0]))),
                             gradient_angle_change_rad=float(wrap(np.arctan2(g1[n,1],g1[n,0])-np.arctan2(g0[n,1],g0[n,0]))))
                    angular=lambda vec,delta:float((vec[0]*delta[1]-vec[1]*delta[0])/max(np.dot(vec,vec),1e-20))
                    oldrow=existing_budget[(i,n)]
                    total_omega=0.
                    for key in oldrow:
                        if key.endswith('_xi') and key not in ('delta_xi',):
                            prefix=key[:-3]; delta=np.array([float(oldrow[key]),float(oldrow[prefix+'_eta'])])
                            value=angular(omid,delta); row['omega_'+prefix+'_rad']=value; total_omega+=value
                            row['omega_'+prefix+'_delta_xi']=float(delta[0]); row['omega_'+prefix+'_delta_eta']=float(delta[1])
                    row['omega_linearization_residual_rad']=row['omega_angle_change_rad']-total_omega
                    total_gradient=0.
                    for name,delta in contributions.items():
                        value=angular(gmid,delta[n]); row['gradient_'+name+'_rad']=value; total_gradient+=value
                        row['gradient_'+name+'_delta_wx']=float(delta[n,0]); row['gradient_'+name+'_delta_wy']=float(delta[n,1])
                    row['gradient_linearization_residual_rad']=row['gradient_angle_change_rad']-total_gradient
                    rates.append(row)
            if i in map_ids:
                saved_maps[i]=dict(fields={name:fields[name] for name in ('w','T','xi','eta','wx','wy','dtv_base','B','Bx','By')},
                                   points=points,coords=coords,time=time,track=tracks[i])
            previous=(i,fields,points)
            print(f'geometry t={time:.1f} ({i})',flush=True)
        metadata=dict(sequence=str(sequence.resolve()),sequence_bytes=sequence.stat().st_size,
                      solver_source_hashes=json.loads(f.attrs['metadata'])['source_sha256'],
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      previous_analysis_source=json.loads((prior/'summary.json').read_text(encoding='utf-8'))['analysis_source_sha256'],
                      phase_windows_s={'before':[2400,2700],'peak':[2700,2880],'after':[2880,3301]},
                      no_new_simulation=True,parcel_ids=list(range(27)))
    write_csv(out/'parcel_geometry.csv',rows); write_csv(out/'orientation_budgets.csv',rates); write_csv(out/'cold_pool.csv',surfaces)
    return rows,rates,surfaces,saved_maps,metadata


def summarize(rows,rates,out):
    trajectories={n:sorted([r for r in rows if r['parcel']==n],key=lambda r:r['time_s']) for n in range(27)}
    nearest=lambda rs,t:min(rs,key=lambda r:abs(r['time_s']-t))
    crossings=[]; decompositions=[]; phases=[]
    for n,rs in trajectories.items():
        for estimator in ('T_product','T','T_trilinear_derivative'):
            # Start from the last positive sample at/before the Eulerian peak.
            positives=[j for j,r in enumerate(rs) if 2700<=r['time_s']<=2800 and r[estimator]>0]
            start=max(positives) if positives else 0
            candidates=[j for j in range(start+1,len(rs)-1) if rs[j-1][estimator]>0 and rs[j][estimator]<0 and rs[j+1][estimator]<0]
            if candidates:
                j=candidates[0]; a,b=rs[j-1:j+1]
                crossings.append(dict(parcel=n,estimator=estimator,t_positive_s=a['time_s'],t_negative_s=b['time_s'],
                                      z_m=b['z_m'],w_ms=b['w'],B_ms2=b['B'],Bthermal_ms2=b['Bthermal'],
                                      dtv_K=b['dtv_base'],surface_dtv_K=b['surface_dtv_K'],
                                      cold_distance_m=b['cold_boundary_distance_m'],angle_deg=b['angle_deg']))
        for name,t0,t1 in (('peak_to_3000',2790,3000),('2700_to_3000',2700,3000),('2790_to_2940',2790,2940)):
            a,b=nearest(rs,t0),nearest(rs,t1)
            result=angular_decomposition([a['xi'],a['eta']],[b['xi'],b['eta']],[a['wx'],a['wy']],[b['wx'],b['wy']])
            if result: decompositions.append(dict(parcel=n,window=name,t0_s=a['time_s'],t1_s=b['time_s'],**result))
        for phase,lo,hi in (('before',2400,2700),('peak',2700,2880),('after',2880,3301),('reversal',2790.25349,2940.14588),('crossing',2790.25349,2820.42587)):
            selected=[r for r in rs if lo<=r['time_s']<hi]
            rr=[r for r in rates if r['parcel']==n and lo<=r['time_s']-r['dt_s']/2<hi]
            phases.append(dict(parcel=n,phase=phase,**{key:float(np.mean([r[key] for r in selected])) for key in
                              ('xi','eta','wx','wy','T_product','Tx_product','Ty_product','angle_deg','w','z_m','B','Bthermal','Bloading','surface_dtv_K','cold_boundary_distance_m')},
                               **{key:float(sum(r[key] for r in rr)) for key in rates[0] if key.endswith(('_rad','_delta_wx','_delta_wy','_delta_xi','_delta_eta'))}))
    write_csv(out/'crossings.csv',crossings); write_csv(out/'geometric_attribution.csv',decompositions); write_csv(out/'phase_sources.csv',phases)
    selected=[r for r in rows if 2400<=r['time_s']<=3301]
    summary=dict(central_samples=[nearest(trajectories[13],t) for t in (2400,2490,2700,2790,2850,2880,2910,2940,3000,3090,3300)],
                 central_crossings=[r for r in crossings if r['parcel']==13],
                 central_attribution=[r for r in decompositions if r['parcel']==13],
                 central_phases=[r for r in phases if r['parcel']==13],
                 crossings={},estimator_sensitivity={})
    for est in ('T_product','T','T_trilinear_derivative'):
        cc=[r for r in crossings if r['estimator']==est]
        summary['crossings'][est]=dict(count=len(cc),first_negative_range_s=[min(r['t_negative_s'] for r in cc),max(r['t_negative_s'] for r in cc)] if cc else [])
        summary['estimator_sensitivity'][est]=dict(sign_agreement_with_product=float(np.mean([np.sign(r[est])==np.sign(r['T_product']) for r in selected])),
                                                  rms_difference=float(np.sqrt(np.mean([(r[est]-r['T_product'])**2 for r in selected]))))
    summary['ensemble_attribution']={name:{key:float(np.median([r[key] for r in decompositions if r['window']==name])) for key in
       ('omega_orientation','gradient_orientation','omega_magnitude','gradient_magnitude','delta_T','delta_alpha_deg','delta_beta_deg','cos_if_only_omega_rotated','cos_if_only_gradient_rotated')}
       for name in ('peak_to_3000','2700_to_3000','2790_to_2940')}
    summary['attribution_counts']={name:dict(
        gradient_orientation_larger=int(sum(abs(r['gradient_orientation'])>abs(r['omega_orientation']) for r in decompositions if r['window']==name)),
        gradient_only_cos_negative=int(sum(r['cos_if_only_gradient_rotated']<0 for r in decompositions if r['window']==name)),
        omega_only_cos_positive=int(sum(r['cos_if_only_omega_rotated']>0 for r in decompositions if r['window']==name)))
        for name in summary['ensemble_attribution']}
    return summary,trajectories,decompositions,phases


def figures(trajectories,decompositions,phases,maps,surfaces,out):
    plt.rcParams.update({'font.size':10,'figure.dpi':140})
    rs=trajectories[13]; t=np.array([r['time_s'] for r in rs])
    fig,axs=plt.subplots(4,2,figsize=(13,13),sharex=True,constrained_layout=True)
    specs=[(('xi','eta'),'Vorticidade horizontal (s⁻¹)'),(('wx','wy'),'Gradiente horizontal de w (s⁻¹)'),
           (('Tx_product','Ty_product','T_product'),'Tilting (s⁻²)'),(('angle_deg',),'Ângulo entre vetores (graus)'),
           (('Bthermal','Bmoist','Bloading','B'),'Flutuabilidade (m s⁻²)'),(('w',),'w (m s⁻¹)'),
           (('surface_dtv_K','dtv_base'),'Anomalia θv: superfície projetada / parcela (K)'),(('z_m',),'Altura da parcela (m)')]
    for ax,(keys,label) in zip(axs.flat,specs):
        for key in keys: ax.plot(t,[r[key] for r in rs],label=key)
        ax.axhline(90 if keys==('angle_deg',) else 0,color='gray',lw=.6); ax.axvline(2790,color='k',ls=':',lw=1)
        ax.set_ylabel(label); ax.grid(alpha=.2); ax.legend(fontsize=8)
    for ax in axs[-1]: ax.set_xlabel('Tempo (s)')
    fig.suptitle('Mesma parcela 13 — geometria, flutuabilidade e trajetória'); fig.savefig(out/'central_geometry.png'); plt.close(fig)
    fig,axs=plt.subplots(3,1,figsize=(11,10),sharex=True,constrained_layout=True)
    for n,rr in trajectories.items():
        for ax,key in zip(axs,('angle_deg','T_product','z_m')): ax.plot(t,[r[key] for r in rr],alpha=.25,color='tab:blue')
    for ax,key in zip(axs,('angle_deg','T_product','z_m')):
        ax.plot(t,[r[key] for r in rs],color='k',label='Parcela 13'); ax.axvline(2790,color='gray',ls=':'); ax.set_ylabel(key); ax.grid(alpha=.2)
    axs[0].axhline(90,color='red',ls='--'); axs[1].axhline(0,color='red',ls='--'); axs[-1].set_xlabel('Tempo (s)')
    fig.suptitle('27 IDs preservados — cruzamento do ângulo de 90°'); fig.savefig(out/'ensemble_geometry.png'); plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for key in ('alpha_deg','beta_deg'): axs[0].plot(t,np.degrees(np.unwrap(np.radians([r[key] for r in rs]))),label=key)
    axs[0].set(xlabel='Tempo (s)',ylabel='Azimute contínuo (graus)',title='Parcela 13: orientação dos dois vetores'); axs[0].legend()
    dec=[r for r in decompositions if r['window']=='peak_to_3000']; keys=['omega_orientation','gradient_orientation','omega_magnitude','gradient_magnitude']
    axs[1].boxplot([[r[k]*1e6 for r in dec] for k in keys],tick_labels=['Orientação ωh','Orientação ∇w','Módulo ωh','Módulo ∇w'])
    axs[1].tick_params(axis='x',rotation=25); axs[1].axhline(0,color='gray'); axs[1].set(ylabel='Contribuição a ΔT (10⁻⁶ s⁻²)',title='2790 → 3000 s: decomposição exata, 27 parcelas')
    fig.savefig(out/'orientation_attribution.png'); plt.close(fig)
    chosen=list(maps.values())[:4]
    fig,axs=plt.subplots(2,len(chosen),figsize=(5*len(chosen),9),constrained_layout=True)
    for col,m in enumerate(chosen):
        x,y,z=m['coords']; p=m['points'][13]; k=int(np.argmin(abs(z-p[2]))); fields=m['fields']
        ix=np.flatnonzero(abs(x-p[0])<4800); iy=np.flatnonzero(abs(y-p[1])<4800)
        xx,yy=np.meshgrid(x[ix]/1000,y[iy]/1000,indexing='ij'); take=lambda key:fields[key][np.ix_(ix,iy,[k])][:,:,0]
        ax=axs[0,col]; im=ax.pcolormesh(xx,yy,take('T')*1e6,cmap='RdBu_r',vmin=-80,vmax=80,shading='auto')
        ww=take('w'); levels=[lev for lev in (-4,-2,0,2,4,8) if ww.min()<lev<ww.max()]
        if levels: ax.contour(xx,yy,ww,levels=levels,colors='gray',linewidths=.7)
        for a,b,color in (('xi','eta','black'),('wx','wy','limegreen')):
            aa,bb=take(a),take(b); norm=np.maximum(np.hypot(aa,bb),1e-12)
            ax.quiver(xx[::2,::2],yy[::2,::2],(aa/norm)[::2,::2],(bb/norm)[::2,::2],color=color,scale=16,width=.005)
        ax.scatter(m['points'][:,0]/1000,m['points'][:,1]/1000,s=7,c='gold',edgecolor='k',linewidth=.2)
        ax.set_title(f"t={m['time']:.0f} s; z={z[k]:.0f} m"); fig.colorbar(im,ax=ax,label='T (10⁻⁶ s⁻²)')
        ax=axs[1,col]; im=ax.pcolormesh(xx,yy,take('B'),cmap='RdBu_r',vmin=-.12,vmax=.12,shading='auto')
        tv=take('dtv_base')
        if tv.min()<-1<tv.max(): ax.contour(xx,yy,tv,levels=[-1],colors='k',linewidths=1.2)
        ax.scatter(p[0]/1000,p[1]/1000,c='gold',edgecolor='k',s=40); fig.colorbar(im,ax=ax,label='B (m s⁻²)')
        for ax in axs[:,col]: ax.set(xlabel='x (km)',ylabel='y (km)',aspect='equal')
    fig.suptitle('Planos na altura da parcela 13: setas unitárias ωh (preto), ∇h w (verde); contornos w em cima, θv′=−1 K em baixo')
    fig.savefig(out/'horizontal_geometry.png'); plt.close(fig)
    fig=plt.figure(figsize=(12,8)); ax=fig.add_subplot(111,projection='3d')
    for n,rr in trajectories.items():
        chosenr=[r for r in rr if 2400<=r['time_s']<=3301]
        ax.plot([r['x_m']/1000 for r in chosenr],[r['y_m']/1000 for r in chosenr],[r['z_m']/1000 for r in chosenr],color='gray',alpha=.25)
    c=ax.scatter([r['x_m']/1000 for r in rs],[r['y_m']/1000 for r in rs],[r['z_m']/1000 for r in rs],c=[r['angle_deg'] for r in rs],cmap='coolwarm',vmin=0,vmax=180,s=25)
    for m in maps.values():
        p=m['points'][13]/1000; r=min(rs,key=lambda r:abs(r['time_s']-m['time']))
        for a,b,color in (('xi','eta','k'),('wx','wy','green')):
            norm=np.hypot(r[a],r[b]); ax.quiver(*p,r[a]/norm,r[b]/norm,0,length=.8,color=color)
        ax.text(*p,f" {m['time']:.0f}s",fontsize=8)
    ax.set(xlabel='x (km)',ylabel='y (km)',zlabel='z (km)',title='Trajetórias preservadas; setas só indicam direção (comprimento arbitrário)')
    fig.colorbar(c,ax=ax,pad=.12,label='Ângulo ωh–∇h w (graus)'); fig.savefig(out/'trajectories_3d.png'); plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(14,6),constrained_layout=True)
    central=[r for r in phases if r['parcel']==13 and r['phase'] in ('crossing','reversal')]
    for ax,prefix in zip(axs,('omega_','gradient_')):
        keys=[k for k in central[0] if k.startswith(prefix) and k.endswith('_rad') and 'angle_change' not in k]
        keys=sorted(keys,key=lambda k:max(abs(r[k]) for r in central),reverse=True)[:7]
        for j,r in enumerate(central): ax.bar(np.arange(len(keys))+(j-1)*.25,[np.degrees(r[k]) for k in keys],width=.25,label=r['phase'])
        ax.set_xticks(np.arange(len(keys)),[k.removeprefix(prefix).removesuffix('_rad') for k in keys],rotation=60,ha='right')
        ax.set(ylabel='Contribuição angular acumulada (graus)',title=prefix+' — parcela 13'); ax.legend(); ax.axhline(0,color='k',lw=.5)
    fig.suptitle('Orçamentos angulares linearizados: crossing 2790–2820; reversal 2790–2940 s')
    fig.savefig(out/'angular_sources.png'); plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for key in ('cold_area_km2','cold_mean_reference_area_km2'):
        axs[0].plot([r['time_s'] for r in surfaces],[r[key] for r in surfaces],label=key)
    axs[0].set_ylabel('Área θv′ < −1 K (km²)'); axs[0].legend()
    axs[1].plot(t,[r['cold_boundary_distance_m']/1000 for r in rs],label='Distância assinada da projeção ao cold pool')
    axs[1].axhline(0,color='k'); axs[1].axvline(2790,color='gray',ls=':'); axs[1].set(xlabel='Tempo (s)',ylabel='Distância (km); negativo = dentro'); axs[1].legend()
    fig.suptitle('Cold pool no nível mais próximo de 100 m; distância de resolução horizontal 600 m')
    fig.savefig(out/'cold_pool_evolution.png'); plt.close(fig)
    fig,axs=plt.subplots(1,len(chosen),figsize=(5*len(chosen),5),constrained_layout=True)
    for ax,m in zip(axs,chosen):
        x,y,z=m['coords']; p=m['points'][13]; j=int(np.argmin(abs(y-p[1]))); ix=np.flatnonzero(abs(x-p[0])<4800); iz=np.flatnonzero(z<=2000)
        xx,zz=np.meshgrid(x[ix]/1000,z[iz]/1000,indexing='ij')
        take=lambda key:m['fields'][key][np.ix_(ix,[j],iz)][:,0,:]
        im=ax.pcolormesh(xx,zz,take('B'),cmap='RdBu_r',vmin=-.12,vmax=.12,shading='auto'); ww=take('w')
        levels=[lev for lev in (-4,-2,0,2,4,8) if ww.min()<lev<ww.max()]
        if levels: ax.contour(xx,zz,ww,levels=levels,colors='k',linewidths=.8)
        ax.scatter(p[0]/1000,p[2]/1000,c='gold',edgecolor='k',s=45)
        ax.set(xlabel='x (km)',ylabel='z (km)',title=f"{m['time']:.0f} s; y={y[j]/1000:.1f} km"); fig.colorbar(im,ax=ax,label='B (m s⁻²)')
    fig.suptitle('Cortes verticais até 2 km: B e contornos w; parcela projetada no plano y mais próximo')
    fig.savefig(out/'vertical_sections.png'); plt.close(fig)


def data_report(summary,trajectories,prior,out):
    track=[r for r in read_csv(prior/'vortex_track.csv') if 2400<=float(r['time_s'])<=3301]
    write_csv(out/'same_component_track.csv',track)
    fig,axs=plt.subplots(2,1,figsize=(11,7),sharex=True,constrained_layout=True)
    rs=trajectories[13]
    for ax,key,trackkey in ((axs[0],'T_product','tilting_at_peak_s2'),(axs[1],'zeta','zeta_low_max_s')):
        ax.plot([r['time_s'] for r in rs],[r[key] for r in rs],label='Mesma parcela 13')
        ax.plot([float(r['time_s']) for r in track],[float(r[trackkey]) for r in track],label='Máximo baixo no mesmo componente (célula pode mudar)')
        ax.set_ylabel(key+' (SI)'); ax.axvline(2790,color='gray',ls=':'); ax.axhline(0,color='gray',lw=.6); ax.legend(fontsize=8)
    axs[-1].set_xlabel('Tempo (s)'); fig.savefig(out/'parcel_vs_component.png'); plt.close(fig)
    lines=['# Tabelas automáticas — inversão do tilting', '',
           'Entrada: sequência existente, somente leitura. Mesmos IDs e coordenadas de `parcels.csv`; nenhuma nova simulação.', '',
           '## Parcela central 13', '',
           '| t (s) | z (m) | ξ (s⁻¹) | η (s⁻¹) | wx (s⁻¹) | wy (s⁻¹) | Tx (s⁻²) | Ty (s⁻²) | ângulo (°) |',
           '|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary['central_samples']:
        lines.append('| '+' | '.join(f'{r[k]:.6g}' for k in ('time_s','z_m','xi','eta','wx','wy','Tx_product','Ty_product','angle_deg'))+' |')
    lines+=['','## Orçamentos angulares da parcela 13 — graus','',
            'Linearização por intervalo de aproximadamente 30 s. Não são experimentos causais; contribuições podem cancelar.', '',
            '| Termo | 2790–2820 s | 2790–2940 s |','|---|---:|---:|']
    phases={r['phase']:r for r in summary['central_phases']}
    for key in phases['crossing']:
        if key.endswith('_rad'):
            lines.append(f"| {key} | {np.degrees(phases['crossing'][key]):.5f} | {np.degrees(phases['reversal'][key]):.5f} |")
    lines+=['','## Sensibilidade','',json.dumps(summary['estimator_sensitivity'],indent=2),'',
            'A primeira amostra negativa exige duas amostras negativas consecutivas depois da última positiva entre 2700 e 2800 s. Tempos são brackets, não instantes exatos.', '',
            '## Arquivos','',
            '`parcel_geometry.csv`: campos, vetores, produtos, ângulos e associação ao cold pool para cada ID/tempo.',
            '`orientation_budgets.csv`: incrementos vetoriais por operador, transporte e resíduos; ângulos em radianos.',
            '`phase_sources.csv`: médias dos campos e somas das fontes por fase; janelas crossing/reversal são adicionais.',
            '`geometric_attribution.csv`: decomposição simétrica exata de ΔT, incluindo módulos; cossenos contrafactuais são geométricos.',
            '`same_component_track.csv`: reprodução do rastreamento anterior, sem redefinir o componente.',
            '`cold_pool.csv`: área no nível 121,66 m, referências ao estado base e à média horizontal.', '',
            'Após aproximadamente 3120 s, o gradiente pode ficar muito pequeno e sua orientação muda rapidamente. O orçamento angular linearizado na fase after completa é mal condicionado; consulte resíduos e incrementos cartesianos. A atribuição principal usa 2790–2940 s.', '',
            'Conclusões e limites físicos: `docs/TILTING_REVERSAL_FINDINGS.md`.']
    (out/'RELATORIO_DADOS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('sequence',type=Path)
    parser.add_argument('--prior',type=Path,required=True); parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    rows,rates,surfaces,maps,metadata=collect(args.sequence,args.prior,args.out)
    summary,trajectories,decompositions,phases=summarize(rows,rates,args.out)
    figures(trajectories,decompositions,phases,maps,surfaces,args.out)
    data_report(summary,trajectories,args.prior,args.out)
    metadata['input_diagnostic_sha256']={name:hashlib.sha256((args.prior/name).read_bytes()).hexdigest() for name in
        ('parcels.csv','lagrangian_budget.csv','vortex_track.csv','summary.json')}
    metadata['current_source_mismatches']=[name for name,digest in metadata['solver_source_hashes'].items()
        if not (ROOT/name).exists() or hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest]
    sha_path=args.sequence.parent/'sequence.sha256.json'
    if sha_path.exists(): metadata['previous_sequence_checksum']=json.loads(sha_path.read_text(encoding='utf-8-sig'))
    metadata['method']=dict(derivatives='NumPy second-order centred differences, nonuniform z; one-sided edge_order=2 at boundaries',
        interpolation='trilinear; same prior parcel positions',g_budget='gradient of native stage delta-w + sampled material transport + explicit residual',
        angles='atan2 orientation; acos relative angle; midpoint angular-source linearization',
        phase_interval_assignment='interval midpoint; no time clipping; end samples for field means',
        cold_pool='theta_v anomaly below -1 K at nearest 100 m; signed grid distance is a horizontal projection, not 3D membership',
        pressure='projection operator retained; no pressure/thermodynamic coupling changes')
    (args.out/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=True),encoding='utf-8')
    (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    print(f'Complete: {args.out}',flush=True)


if __name__=='__main__': main()
