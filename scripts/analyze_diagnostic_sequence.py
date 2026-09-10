"""Independent Eulerian audit, overlap-based vortex identity, and RK4 parcels.

All plots and the numerical report are regenerated from sequence.h5. Never
imports solver differentiation, vorticity budget, or classification routines.
The recorder's standalone NumPy derivative functions are reused for consistent
linear curl of recorded increments; they are independent of solver stencils.
"""
from pathlib import Path
import argparse
import csv
import json
import hashlib
import zipfile
import sys
import numpy as np
import h5py
from scipy.ndimage import label
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from storm_dynamics.diagnostic_capture import curl, centered, kinematics


def write_csv(path,rows):
    if not rows: return
    keys=list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=keys); w.writeheader(); w.writerows(rows)


def native(group,nk=None):
    return tuple(group[n][...] if nk is None else group[n][:,:,:nk+(n=='w')] for n in ('u','v','w'))


def candidate_mask(zeta,z,threshold=.003):
    mask=(zeta>=threshold)&(z[None,None,:]<=2000.)
    mask[:2]=False; mask[-2:]=False; mask[:,:2]=False; mask[:,-2:]=False
    return label(mask)  # 6-neighbour spatial continuity, no intensity ratio proxy


def associate(labels,old,old_center,coordinates,weights):
    """Choose only overlapping components; never jump to a new domain maximum."""
    ids,counts=np.unique(labels[old],return_counts=True)
    valid=[(int(n),int(c)) for n,c in zip(ids,counts) if n and c>=3]
    if not valid: return None,'lost',0.
    ranked=[]
    for n,c in valid:
        new=labels==n
        score=c/float(np.count_nonzero(old|new))
        ranked.append((score,n))
    ranked.sort(reverse=True)
    score,n=ranked[0]
    if score<.02: return None,'lost_low_overlap',score
    status='ambiguous' if len(ranked)>1 and ranked[1][0]>=.8*score else 'overlap'
    if status=='ambiguous': return None,status,score
    return labels==n,status,score


def audit_and_track(f,out,seed_time=1800.):
    groups=list(f['snapshots'].values()); times=np.array([g.attrs['time_s'] for g in groups])
    nk=int(f.attrs['budget_nz_with_halo'])
    coords=tuple(f['grid'][n][:] for n in ('xc','yc','zc')); lowcoords=(*coords[:2],coords[2][:nk])
    z=lowcoords[2]; dz=np.diff(f['grid/zf'][:])[:nk]
    volume=np.broadcast_to(dz[None,None,:],(len(coords[0]),len(coords[1]),nk)).copy()
    interior=np.zeros(volume.shape,bool); interior[2:-2,2:-2,:]=(z<=2000)[None,None,:]
    row_budget=[]; zetas=[]; ws=[]; labels_all=[]; previous=None
    pressure_ranges=[]; operator_max={}
    for i,g in enumerate(groups):
        vel=native(g,nk); kin=kinematics(vel,lowcoords); om=kin['omega']
        zetas.append(om[2]); ws.append(centered(vel)[2])
        labels_all.append(candidate_mask(om[2],z)[0])
        pressure_ranges.append(dict(time_s=float(times[i]),p_min=float(g['p'][:].min()),p_max=float(g['p'][:].max()),
                                    p_dyn_min=float(g['p_dyn'][:].min()) if 'p_dyn' in g else None,
                                    p_dyn_max=float(g['p_dyn'][:].max()) if 'p_dyn' in g else None))
        if previous is not None:
            duration=times[i]-times[i-1]
            increments={}
            for name,stage in g['increments'].items():
                changes=native(stage)
                operator_max[name]=max(operator_max.get(name,0.),max(float(abs(a).max()) for a in changes))
                increments[name]=curl(changes,lowcoords)
            total=sum(increments.values()); delta=om-previous
            error=delta-total
            integral=g['kinematic_integrals']
            physical=sum(integral[n][:] for n in ('vector_advection','vector_stretching','vector_dilatation'))
            mismatch=increments['advection']-physical
            rms=lambda arr:float(np.sqrt(np.average(arr[interior]**2,weights=volume[interior])))
            u,v,w=vel; rho=f['base/rho0'][:nk]; rhow=f['base/rho0_wface'][:nk+1]
            divmass=(np.diff(u,axis=0)*rho/np.diff(coords[0]).mean()+np.diff(v,axis=1)*rho/np.diff(coords[1]).mean()+
                     np.diff(w*rhow,axis=2)/dz)
            row_budget.append(dict(time_s=times[i],dt_interval_s=duration,
                                   closure_max_s2=float(np.abs(error[:,interior]).max()/duration),
                                   closure_zeta_rms_s2=rms(error[2])/duration,
                                   observed_zeta_rms_s2=rms(delta[2])/duration,
                                   continuum_advection_mismatch_rms_s2=rms(mismatch[2])/duration,
                                   advective_curl_rms_s2=rms(increments['advection'][2])/duration))
            row_budget[-1]['native_mass_divergence_max_kg_m3_s']=float(abs(divmass[interior]).max())
        previous=om
        if i%20==0: print(f'Eulerian {i}/{len(groups)} t={times[i]:.1f}',flush=True)
    anchor=int(np.argmin(abs(times-seed_time)))
    region=(z>=500)&(z<=1500)
    score=zetas[anchor]*np.maximum(ws[anchor],0)
    score=np.where(interior & region[None,None,:] & (labels_all[anchor]>0),score,-np.inf)
    ijk=np.unravel_index(np.argmax(score),score.shape)
    if not np.isfinite(score[ijk]): raise RuntimeError('No rotating updraft candidate at seed time; no identity fabricated')
    masks={anchor:labels_all[anchor]==labels_all[anchor][ijk]}; statuses={anchor:('seed',1.)}
    for direction in (-1,1):
        old=masks[anchor]
        for i in range(anchor+direction,len(groups) if direction>0 else -1,direction):
            mask,status,overlap=associate(labels_all[i],old,None,lowcoords,None)
            statuses[i]=(status,overlap)
            if mask is None: break
            masks[i]=mask; old=mask
    rows=[]; profiles=[]; maps={}; xx,yy=np.meshgrid(coords[0],coords[1],indexing='ij')
    for i,g in enumerate(groups):
        if i not in masks:
            rows.append(dict(index=i,time_s=times[i],status=statuses.get(i,('untracked',0))[0])); continue
        mask=masks[i]; zz=zetas[i]; ww=ws[i]; weight=zz*volume*mask; mass=weight.sum()
        local_kinematics=kinematics(native(g,nk),lowcoords)
        cx=float((weight*xx[:,:,None]).sum()/mass); cy=float((weight*yy[:,:,None]).sum()/mass)
        iz=np.where(mask.any(axis=(0,1)))[0]; low=mask & (z<=500)[None,None,:]
        target=low if low.any() else mask
        pos=np.unravel_index(np.argmax(np.where(target,zz,-np.inf)),zz.shape)
        r=np.hypot(xx-coords[0][pos[0]],yy-coords[1][pos[1]])
        ring=(r>=3000)&(r<=6000); core=r<=1200
        p=g['p_dyn'][:,:,pos[2]] if 'p_dyn' in g else g['p'][:,:,pos[2]]
        pressure=float(np.median(p[core])-np.median(p[ring])) if ring.any() else np.nan
        cloud=g['ql'][:,:,:nk]+g['qi'][:,:,:nk]
        cindices=np.where((cloud>1e-5).any(axis=(0,1)))[0]
        inside=(cloud>1e-5)&mask
        baseidx=np.where(inside.any(axis=(0,1)))[0]
        row=dict(index=i,time_s=times[i],status=statuses[i][0],overlap_iou=statuses[i][1],
                 center_x_m=cx,center_y_m=cy,zeta_low_max_s=float(zz[pos]),peak_z_m=float(z[pos[2]]),
                 w_at_peak_ms=float(ww[pos]),zeta_volume_mean_s=float(np.sum(zz*volume*mask)/np.sum(volume*mask)),
                 stretching_at_peak_s2=float(local_kinematics['stretching'][pos]),
                 tilting_at_peak_s2=float(local_kinematics['tilting'][pos]),
                 convergence_at_peak_s=float(local_kinematics['convergence'][pos]),
                 dilatation_at_peak_s2=float(local_kinematics['vector_dilatation'][2][pos]),
                 z_min_m=float(z[iz[0]]),z_max_m=float(z[iz[-1]]),
                 touches_lowest_resolved_level=bool(mask[:,:,0].any()),
                 core_pressure_minus_ring_Pa=pressure,
                 cloud_base_in_component_m=float(z[baseidx[0]]) if len(baseidx) else np.nan,
                 domain_cloud_base_m=float(z[cindices[0]]) if len(cindices) else np.nan,
                 lowest_level_cloud_in_component=bool(inside[:,:,0].any()),
                 peak_i=int(pos[0]),peak_j=int(pos[1]),peak_k=int(pos[2]))
        rows.append(row)
        for k in iz:
            weight2=zz[:,:,k]*mask[:,:,k]
            profiles.append(dict(time_s=times[i],z_m=float(z[k]),
                                 x_m=float((weight2*xx).sum()/weight2.sum()),y_m=float((weight2*yy).sum()/weight2.sum()),
                                 zeta_max_s=float(zz[:,:,k][mask[:,:,k]].max()),
                                 w_mean_ms=float(np.average(ww[:,:,k][mask[:,:,k]],weights=weight2[mask[:,:,k]]))))
    write_csv(out/'eulerian_budget.csv',row_budget); write_csv(out/'vortex_track.csv',rows)
    write_csv(out/'pressure_ranges.csv',pressure_ranges)
    write_csv(out/'vertical_coherence.csv',profiles)
    valid=[r for r in rows if r.get('peak_z_m',np.inf)<=500 and r['status']!='ambiguous']
    if not valid: raise RuntimeError('Tracked component never reaches <=500m; no low-level parcel seed fabricated')
    peak=max(valid,key=lambda r:r['zeta_low_max_s']); peak_index=peak['index']
    # 27 deterministic seeds around the tracked low-level maximum; remove
    # points outside domain, retain weak/descending parcels (no success filter).
    point=np.array([coords[0][peak['peak_i']],coords[1][peak['peak_j']],coords[2][peak['peak_k']]])
    offsets=np.array(np.meshgrid([-300.,0,300.],[-300.,0,300.],[-30.,0,30.],indexing='ij')).reshape(3,-1).T
    seeds=point+offsets
    seeds=seeds[(seeds[:,2]>=coords[2][0]) & (seeds[:,2]<=2000)]
    # Memory-heavy masks/labels are not held during trajectory interpolation.
    summary=dict(seed_time_s=float(times[anchor]),seed_location_m=[float(lowcoords[a][ijk[a]]) for a in range(3)],
                 threshold_s=.003,component_connectivity=6,seed_selection='max positive zeta*w, 500–1500 m at requested 1800 s',
                 identity_rule='spatial component overlap only; IoU >= 0.02 and >=3 cells; competing overlap >=80% is ambiguous',
                 tracked_start_s=float(min(times[i] for i in masks)),tracked_end_s=float(max(times[i] for i in masks)),
                 peak=peak,closure_max_s2=max(r['closure_max_s2'] for r in row_budget),
                 closure_relative_max=max(r['closure_zeta_rms_s2']/max(r['observed_zeta_rms_s2'],1e-30) for r in row_budget),
                 continuum_mismatch_relative_median=float(np.median([r['continuum_advection_mismatch_rms_s2']/max(r['advective_curl_rms_s2'],1e-30) for r in row_budget])))
    summary['mass_divergence_max_kg_m3_s']=max(r['native_mass_divergence_max_kg_m3_s'] for r in row_budget)
    summary['identity_stop_events']=[dict(time_s=float(times[i]),status=status,overlap_iou=float(overlap))
                                      for i,(status,overlap) in statuses.items() if i not in masks]
    summary['operator_interval_increment_max_ms']=operator_max
    summary['thermodynamic_p_identically_zero']=all(r['p_min']==0 and r['p_max']==0 for r in pressure_ranges)
    after=[r for r in valid if r['index']>peak_index]
    weakening=None
    for a,b,c in zip(after,after[1:],after[2:]):
        if b['index']==a['index']+1 and c['index']==b['index']+1 and max(r['zeta_low_max_s'] for r in (a,b,c))<.8*peak['zeta_low_max_s']:
            weakening=a; break
    summary['sustained_20percent_weakening']=weakening
    summary['visible_lowest_level_cloud_frames']=sum(r.get('lowest_level_cloud_in_component',False) for r in rows)
    # Separate cloud condensate, velocity-derived rotation, and dynamic pressure
    # in a vertical section through the *same* tracked peak.
    group=groups[peak_index]; j=peak['peak_j']; x=coords[0]/1000; heights=z/1000
    cloud=(group['ql'][:,j,:nk]+group['qi'][:,j,:nk])*1000
    pressure=group['p_dyn'][:,j,:nk] if 'p_dyn' in group else group['p'][:,j,:nk]
    pressure=pressure-np.median(pressure,axis=0,keepdims=True)
    fig,axes=plt.subplots(1,3,figsize=(13,5),constrained_layout=True)
    for ax,data,title,unit in zip(axes,(zetas[peak_index][:,j,:],cloud,pressure),('Reconstructed zeta','Cloud ql + qi','Dynamic pressure − section median'),('s$^{-1}$','g/kg','Pa')):
        im=ax.pcolormesh(x,heights,data.T,shading='auto',cmap='RdBu_r' if title!='Cloud ql + qi' else 'Blues')
        ax.set(xlim=(point[0]/1000-10,point[0]/1000+10),ylim=(0,2),xlabel='x (km)',ylabel='z (km)',title=title)
        fig.colorbar(im,ax=ax,label=unit)
    fig.suptitle(f'Same section at t={times[peak_index]:.1f} s; cloud is not a pressure proxy')
    fig.savefig(out/'pressure_vs_visible_cloud.png',dpi=150); plt.close(fig)
    # Figures: a tracked component, no connection from independent maxima.
    fig,axes=plt.subplots(3,1,figsize=(10,10),sharex=True,constrained_layout=True)
    tr=[r for r in rows if 'zeta_low_max_s' in r]; t=[r['time_s']/60 for r in tr]
    plotted=lambda key:[r[key] if r['peak_z_m']<=500 else np.nan for r in tr]
    axes[0].plot(t,plotted('zeta_low_max_s')); axes[0].set_ylabel('Tracked peak below 500 m (s$^{-1}$)')
    axes[1].plot(t,plotted('w_at_peak_ms')); axes[1].axhline(0,color='grey'); axes[1].set_ylabel('w at tracked peak (m/s)')
    axes[2].plot(t,plotted('core_pressure_minus_ring_Pa')); axes[2].set_ylabel('Dynamic core − ring (Pa)'); axes[2].set_xlabel('Simulation time (min)')
    for ax in axes: ax.grid(alpha=.3)
    fig.suptitle('One overlap-tracked component: pressure is independent of visibility')
    fig.savefig(out/'tracked_vortex.png',dpi=150); plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4),constrained_layout=True)
    for name in ('stretching_at_peak_s2','tilting_at_peak_s2','dilatation_at_peak_s2'):
        ax.plot(t,plotted(name),label=name)
    ax.axhline(0,color='grey'); ax.set(xlabel='Simulation time (min)',ylabel='Local term at tracked peak (s$^{-2}$)')
    ax.legend(); ax.grid(alpha=.3); fig.savefig(out/'tracked_local_terms.png',dpi=150); plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
    ax.semilogy([r['time_s']/60 for r in row_budget],[max(r['closure_zeta_rms_s2'],1e-25) for r in row_budget],label='Discrete closure residual')
    ax.semilogy([r['time_s']/60 for r in row_budget],[r['continuum_advection_mismatch_rms_s2'] for r in row_budget],label='Actual advection curl − continuous terms')
    ax.set(xlabel='Time (min)',ylabel='RMS tendency (s$^{-2}$)'); ax.legend(); ax.grid(alpha=.3)
    fig.savefig(out/'budget_closure.png',dpi=150); plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4),constrained_layout=True)
    scatter=ax.scatter([p['time_s']/60 for p in profiles],[p['z_m']/1000 for p in profiles],c=[p['zeta_max_s'] for p in profiles],s=12)
    ax.set(xlabel='Time (min)',ylabel='Height (km)',ylim=(0,2)); fig.colorbar(scatter,ax=ax,label='zeta in same connected component (s$^{-1}$)')
    fig.savefig(out/'vertical_coherence.png',dpi=150); plt.close(fig)
    return summary,peak_index,seeds,rows


class VelocityInterval:
    def __init__(self,left,right,f):
        self.t0=float(left.attrs['time_s']); self.t1=float(right.attrs['time_s'])
        coords={n:f['grid'][n][:] for n in ('xc','yc','zc','xf','yf','zf')}
        # Compact causal-restart files retain only the diagnosed lower column.
        # Match coordinate vectors to the native staggered arrays actually saved.
        coords['zc']=coords['zc'][:left['u'].shape[2]]
        coords['zf']=coords['zf'][:left['w'].shape[2]]
        axes=((coords['xf'],coords['yc'],coords['zc']), (coords['xc'],coords['yf'],coords['zc']), (coords['xc'],coords['yc'],coords['zf']))
        self.interpolators=[[RegularGridInterpolator(ax,g[n][...],bounds_error=False,fill_value=np.nan)
                             for ax,n in zip(axes,('u','v','w'))] for g in (left,right)]
        self.low=np.array([coords[n][0] for n in ('xc','yc','zc')]); self.high=np.array([coords[n][-1] for n in ('xc','yc','zc')])

    def velocity(self,t,points):
        alpha=(t-self.t0)/(self.t1-self.t0)
        return sum(weight*np.stack([ip(points) for ip in ips],axis=-1) for weight,ips in zip((1-alpha,alpha),self.interpolators))


def rk4(interval,points,start,end,max_step):
    steps=int(np.ceil(abs(end-start)/max_step)); dt=(end-start)/steps; p=points.copy(); t=start
    for _ in range(steps):
        a=interval.velocity(t,p); b=interval.velocity(t+dt/2,p+dt*a/2)
        c=interval.velocity(t+dt/2,p+dt*b/2); d=interval.velocity(t+dt,p+dt*c)
        p+=dt*(a+2*b+2*c+d)/6; t+=dt
        bad=((p<interval.low)|(p>interval.high)).any(axis=1); p[bad]=np.nan
    return p


def trajectories(f,seed_index,seeds,max_step=2.,stride=1):
    groups=list(f['snapshots'].values()); ids=sorted(set(range(0,len(groups),stride))|{seed_index,len(groups)-1})
    positions={seed_index:seeds.copy()}
    for direction in (-1,1):
        start=ids.index(seed_index); current=seeds.copy()
        for j in range(start,len(ids)-1) if direction==1 else range(start,0,-1):
            a,b=(ids[j],ids[j+1]) if direction==1 else (ids[j-1],ids[j])
            interval=VelocityInterval(groups[a],groups[b],f)
            current=rk4(interval,current,interval.t0 if direction==1 else interval.t1,interval.t1 if direction==1 else interval.t0,max_step)
            positions[b if direction==1 else a]=current.copy()
    return positions


def lagrangian_budget(f,positions,out):
    groups=list(f['snapshots'].values()); nk=int(f.attrs['budget_nz_with_halo'])
    coords=tuple(f['grid'][n][:] for n in ('xc','yc','zc')); coords=(*coords[:2],coords[2][:nk])
    sample=lambda a,p:RegularGridInterpolator(coords,a,bounds_error=False,fill_value=np.nan)(p)
    rows=[]; budget=[]; previous=None
    for i,g in enumerate(groups):
        p=positions[i]; om=curl(native(g,nk),coords)
        values=np.stack([sample(a,p) for a in om],axis=-1)
        for n,(point,val) in enumerate(zip(p,values)):
            rows.append(dict(index=i,parcel=n,time_s=float(g.attrs['time_s']),x_m=point[0],y_m=point[1],z_m=point[2],
                             xi_s=val[0],eta_s=val[1],zeta_s=val[2],
                             budget_valid=bool(np.isfinite(val).all() and point[2]<=2000 and point[2]>=coords[2][0])))
        if previous is not None:
            # Mid-path interval quadrature of dt-integrated Eulerian fields.
            # This is NOT a machine-precision material closure claim.
            midpoint=(positions[i]+positions[i-1])/2
            inside=(midpoint[:,2]<=2000)&(positions[i][:,2]<=2000)&(positions[i-1][:,2]<=2000)
            ki=g['kinematic_integrals']
            production={name:np.stack([sample(a,midpoint) for a in ki[name][:]],axis=-1)
                        for name in ('vector_stretching','vector_dilatation')}
            production['axial_stretching']=sample(ki['stretching'][:],midpoint)
            production['tilting']=sample(ki['tilting'][:],midpoint)
            adv=curl(native(g['increments/advection']),coords)
            continuous=sum(ki[name][:] for name in ('vector_advection','vector_stretching','vector_dilatation'))
            production['advection_operator_remainder']=np.stack([sample(a,midpoint) for a in adv-continuous],axis=-1)
            for name,stage in g['increments'].items():
                if name!='advection': production[name]=np.stack([sample(a,midpoint) for a in curl(native(stage),coords)],axis=-1)
            total=sum(v for k,v in production.items() if k not in ('axial_stretching','tilting'))
            change=values-previous
            for n in range(len(p)):
                row=dict(index=i,parcel=n,time_s=float(g.attrs['time_s']),valid=bool(inside[n] and np.isfinite(total[n]).all() and np.isfinite(change[n]).all()))
                for axis,name in enumerate(('xi','eta','zeta')):
                    row['delta_'+name]=change[n,axis]; row['residual_'+name]=change[n,axis]-total[n,axis]
                    for key,val in production.items():
                        if val.ndim==2: row[key+'_'+name]=val[n,axis]
                row['stretching_zeta']=production['axial_stretching'][n]; row['tilting_zeta']=production['tilting'][n]
                budget.append(row)
        previous=values
        if i%20==0: print(f'Lagrangian budget {i}/{len(groups)}',flush=True)
    write_csv(out/'parcels.csv',rows); write_csv(out/'lagrangian_budget.csv',budget)
    valid=[r for r in budget if r['valid']]
    rms=lambda key:float(np.sqrt(np.mean([r[key]**2 for r in valid]))) if valid else None
    fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True)
    for n in sorted({r['parcel'] for r in rows}):
        rr=[r for r in rows if r['parcel']==n]
        axes[0].plot([r['time_s']/60 for r in rr],[r['z_m']/1000 for r in rr],alpha=.6)
        axes[1].plot([r['time_s']/60 for r in rr],[r['zeta_s'] for r in rr],alpha=.6)
    axes[0].set(xlabel='Time (min)',ylabel='Parcel height (km)'); axes[1].set(xlabel='Time (min)',ylabel='Parcel zeta (s$^{-1}$)')
    for ax in axes: ax.grid(alpha=.3)
    fig.savefig(out/'lagrangian_paths.png',dpi=150); plt.close(fig)
    # Central seed's terms, no averaging of distinct histories into a trajectory.
    parcel=len(positions[0])//2; rr=[r for r in budget if r['parcel']==parcel and r['valid']]
    fig,ax=plt.subplots(figsize=(10,5),constrained_layout=True)
    for key in ('delta_zeta','stretching_zeta','tilting_zeta','surface_drag_zeta','les_zeta','projection_zeta','advection_operator_remainder_zeta','residual_zeta'):
        ax.plot([r['time_s']/60 for r in rr],np.cumsum([r[key] for r in rr]),label=key)
    ax.set(xlabel='Time (min); only valid intervals included',ylabel='Cumulative change (s$^{-1}$)'); ax.legend(fontsize=8,ncol=2); ax.grid(alpha=.3)
    fig.savefig(out/'lagrangian_budget.png',dpi=150); plt.close(fig)
    return dict(parcel_count=len(positions[0]),valid_budget_intervals=len(valid),total_budget_intervals=len(budget),
                lagrangian_delta_rms_s=rms('delta_zeta'),lagrangian_residual_rms_s=rms('residual_zeta'),
                records=rows,budget=budget)


def write_report(summary,result,track,out,f):
    """Evidence labels are driven by measured residuals, never desired outcomes."""
    peak=summary['peak']; groups=list(f['snapshots'].values())
    central=result['parcel_count']//2
    parcel=[r for r in result['records'] if r['parcel']==central]
    intervals={r['index']:r for r in result['budget'] if r['parcel']==central}
    seed=peak['index']; start=seed; end=seed
    while start in intervals and intervals[start]['valid']: start-=1
    while end+1 in intervals and intervals[end+1]['valid']: end+=1
    continuous=[intervals[i] for i in range(start+1,end+1)]
    origin=parcel[start]; last=parcel[end]
    terms=('stretching_zeta','tilting_zeta','vector_dilatation_zeta','surface_drag_zeta','les_zeta','projection_zeta',
           'coriolis_zeta','buoyancy_zeta','advection_operator_remainder_zeta','residual_zeta')
    totals={k:float(sum(r.get(k,0) for r in continuous)) for k in terms}
    before=[r for r in continuous if r['index']<=seed]
    after=[r for r in continuous if r['index']>seed]
    sums=lambda rr:{k:float(sum(r.get(k,0) for r in rr)) for k in terms}
    summary['representative_parcel']=dict(id=central,origin=origin,last=last,contiguous_interval_count=len(continuous),
                                         all_contributions_s=totals,before_vortex_peak=sums(before),after_vortex_peak=sums(after))
    fig,axes=plt.subplots(1,3,figsize=(15,5),constrained_layout=True)
    for ax,component in zip(axes,('xi','eta','zeta')):
        for name in ('delta','vector_stretching','vector_dilatation','buoyancy','surface_drag','les','projection','advection_operator_remainder','residual'):
            ax.plot([r['time_s']/60 for r in continuous],np.cumsum([r.get(name+'_'+component,0.) for r in continuous]),label=name)
        ax.set(xlabel='Time (min)',ylabel='Cumulative '+component+' change (s$^{-1}$)'); ax.grid(alpha=.3)
    axes[-1].legend(fontsize=7); fig.suptitle('Representative parcel: one contiguous 0–2 km segment')
    fig.savefig(out/'horizontal_vorticity_sources.png',dpi=150); plt.close(fig)
    origins=[]
    parcel_phase_rows=[]
    for n in range(result['parcel_count']):
        rr=[r for r in result['records'] if r['parcel']==n]
        ii={r['index']:r for r in result['budget'] if r['parcel']==n}
        s=seed
        while s in ii and ii[s]['valid']: s-=1
        origins.append(dict(parcel=n,**{k:rr[s][k] for k in ('time_s','x_m','y_m','z_m','xi_s','eta_s','zeta_s')}))
        e=seed
        while e+1 in ii and ii[e+1]['valid']: e+=1
        for phase,indices in (('before_peak',range(s+1,seed+1)),('after_peak',range(seed+1,e+1))):
            selected=[ii[i] for i in indices]
            parcel_phase_rows.append(dict(parcel=n,phase=phase,intervals=len(selected),
                                          start_s=rr[s if phase=='before_peak' else seed]['time_s'],
                                          end_s=rr[seed if phase=='before_peak' else e]['time_s'],
                                          delta_zeta=float(sum(r['delta_zeta'] for r in selected)),**sums(selected)))
    write_csv(out/'parcel_origins_contiguous.csv',origins)
    write_csv(out/'parcel_phase_contributions.csv',parcel_phase_rows)
    summary['parcel_ensemble_phase']={}
    for phase in ('before_peak','after_peak'):
        selected=[r for r in parcel_phase_rows if r['phase']==phase and r['intervals']>0]
        summary['parcel_ensemble_phase'][phase]=dict(count=len(selected),negative_delta_count=sum(r['delta_zeta']<0 for r in selected),
                                                    negative_tilting_count=sum(r['tilting_zeta']<0 for r in selected),
                                                    medians={key:float(np.median([r[key] for r in selected])) if selected else None for key in ('delta_zeta',*terms)})
    bz=f['grid/zc'][:]
    environment_xi=float(np.interp(origin['z_m'],bz,-np.gradient(f['base/v0'][:],bz,edge_order=2))) if 'v0' in f['base'] else np.nan
    environment_eta=float(np.interp(origin['z_m'],bz,np.gradient(f['base/u0'][:],bz,edge_order=2))) if 'u0' in f['base'] else np.nan
    summary['representative_parcel']['base_horizontal_vorticity_at_origin_s']=[environment_xi,environment_eta]
    rms=result['lagrangian_residual_rms_s']; delta=result['lagrangian_delta_rms_s']
    relative=rms/max(delta,1e-30) if rms is not None else None
    summary['lagrangian']['relative_rms_residual']=relative
    weakness=summary['sustained_20percent_weakening']
    weakened=(f"A primeira queda de pelo menos 20% sustentada por três saídas consecutivas começa em {weakness['time_s']:.1f} s "
              f"(zeta={weakness['zeta_low_max_s']:.5f} s^-1)." if weakness else
              'Não foi observada queda de 20% sustentada por três saídas consecutivas durante a identidade rastreada; perda de identidade não é enfraquecimento comprovado.')
    table='\n'.join(f'| {key} | {value:+.6g} |' for key,value in totals.items())
    phase_table='\n'.join(f"| {key} | {sums(before)[key]:+.6g} | {sums(after)[key]:+.6g} |" for key in terms)
    settings=summary['metadata']['config']; simcfg=settings['sim']
    text=f'''# Relatório diagnóstico temporal — sequência instrumentada

Nenhuma equação, parametrização ou termodinâmica foi alterada. A sequência foi integrada desde a condição analítica inicial, sem importar um cache histórico. As conclusões abaixo são deste caso idealizado e desta malha; não constituem reprodução observacional de Moore nem diagnóstico de todos os ninhos antigos.

## Simulação e produtos

- Intervalo: {groups[0].attrs['time_s']:.1f}–{groups[-1].attrs['time_s']:.1f} s; {len(groups)} saídas sincronizadas, com horários reais e todos os dt nativos registrados.
- Malha: {simcfg['grid']['nx']} × {simcfg['grid']['ny']} × {simcfg['grid']['nz']}; dx={np.diff(f['grid/xc'][:]).mean():.1f} m; primeiro centro z={f['grid/zc'][0]:.2f} m. Não é uma malha que resolva um funil estreito.
- Domínio fixo em referencial relativo à tempestade; contornos laterais periódicos; configuração da superfície copiada do caso existente `tornado_intensity_L_gpu.py`, sem ajustes motivados pelo resultado.
- Arquivo `../sequence.h5`: velocidades nativas e estado termodinâmico tridimensionais na coluna completa; pressões `p` e `p_dyn` separadas; estado base e coordenadas; incrementos nativos por operador e integrais cinemáticas nos primeiros 2 km, com três níveis de halo.
- `../metadata.json`: configuração completa, versão do Python/NumPy, equipamento, SHA-256 dos fontes, revisão git e diff executado. O arquivo HDF5 contém também os metadados iniciais.
- CSVs: orçamento Euleriano, rastreamento, coerência vertical, trajetórias, origens delimitadas e orçamento lagrangiano. Figuras PNG são geradas automaticamente pelo analisador.

## O que está comprovado nesta execução

1. **A captura é passiva nos testes realizados.** Testes CPU/GPU comparam bit a bit campos prognósticos; há verificação específica do caminho de baixa memória e de suas duas pressões. Isso verifica a instrumentação nos casos de teste, não é uma segunda integração completa deste caso.
2. **Os incrementos discretos fecham.** Máximo residual do rotacional da mudança menos a soma dos rotacionais dos incrementos: {summary['closure_max_s2']:.3e} s^-2. Maior razão entre RMS residual e RMS da mudança de zeta: {summary['closure_relative_max']:.3e}. Esse fechamento telescópico detecta mudanças de velocidade não contabilizadas; sozinho não valida uma equação contínua.
3. **A continuidade nativa foi medida.** Máximo interior abaixo de 2 km de |div(rho0 u)| nos instantes salvos: {summary['mass_divergence_max_kg_m3_s']:.3e} kg m^-3 s^-1. São diferenças nas faces e na malha original, não divergência de velocidades recentradas.
4. **A pressão dinâmica e a visibilidade são campos distintos.** `p` termodinâmico permaneceu identicamente zero: {summary['thermodynamic_p_identically_zero']}. `p_dyn` foi salvo separadamente. Houve {summary['visible_lowest_level_cloud_frames']} saídas com ql+qi > 1e-5 kg/kg no nível mais baixo dentro do componente rastreado. Ausência de condensado não equivale a ausência de circulação; condensado em um componente amplo tampouco confirma um funil estreito.

## O que é suportado pelo rastreamento

Foi escolhido um componente ciclônico de zeta >= 0.003 s^-1, pelo maior zeta × w positivo entre 500 e 1500 m no instante próximo de {summary['seed_time_s']:.1f} s. A escolha foi predefinida no analisador. Depois da semente, só componentes sobrepostos são aceitos. Não há procura por outro máximo do domínio após a perda. Competidores com escores próximos são marcados ambíguos. A conexão espacial usa seis vizinhos na mesma malha.

A identidade espacial foi acompanhada de {summary['tracked_start_s']:.1f} a {summary['tracked_end_s']:.1f} s. Seu pico abaixo de 500 m ocorreu em **{peak['time_s']:.1f} s**, com zeta={peak['zeta_low_max_s']:.5f} s^-1, z={peak['peak_z_m']:.1f} m e w={peak['w_at_peak_ms']:.3f} m/s. {weakened}

Eventos de interrupção da identidade: {summary['identity_stop_events']}. Nenhum ramo é escolhido depois de uma associação ambígua. As parcelas continuam sendo advectadas individualmente, sem exigir que permaneçam no componente Euleriano.

No pico, a pressão dinâmica média por mediana no disco de 1.2 km, relativa à mediana do anel de 3–6 km na mesma altura, foi {peak['core_pressure_minus_ring_Pa']:+.2f} Pa. É uma anomalia de pressão da tempestade em torno do pico rastreado; não foi isolada a contribuição exclusivamente centrífuga. A figura `pressure_vs_visible_cloud.png` mantém pressão, zeta e condensado em painéis distintos.

O critério de sobreposição é Euleriano e pode seguir uma estrutura ampla em fusões ou cisalhamento. Ele não identifica um tubo material único, e o valor de 0.003 s^-1 é um limiar de análise, não de classificação de tornado. Altura mínima e máxima do mesmo componente estão em `vertical_coherence.csv`; não foram ligadas colunas de máximos independentes.

## Origem e evolução lagrangiana resolvida

Foram lançadas {result['parcel_count']} sementes determinísticas em torno do pico rastreado abaixo de 500 m, com deslocamentos horizontais de ±300 m e verticais de ±30 m. Não se filtraram parcelas descendentes ou de fraco crescimento. A advecção das parcelas usa u/v/w nas faces nativas, interpolação trilinear no espaço, linear entre saídas e RK4 para trás e para frente no tempo.

O orçamento é válido apenas em segmentos contíguos dentro da faixa espacial amostrada e abaixo de 2 km. Parcelas que deixam os centros verticais resolvidos ou os limites laterais deixam de ser extrapoladas. Não se aplicou continuação periódica artificial das trajetórias. A origem tabulada é a primeira posição **resolvida e contígua** até a semente, não uma alegação de origem física antes dessa posição.

Para a parcela representativa {central}, o segmento contíguo é {origin['time_s']:.1f}–{last['time_s']:.1f} s. A origem resolvida fica em (x,y,z)=({origin['x_m']:.1f}, {origin['y_m']:.1f}, {origin['z_m']:.1f}) m, com (xi,eta,zeta)=({origin['xi_s']:.4g}, {origin['eta_s']:.4g}, {origin['zeta_s']:.4g}) s^-1. As demais origens, sem seleção posterior por sucesso, estão em `parcel_origins_contiguous.csv`.

As mudanças acumuladas de zeta nesse segmento são:

| Termo | Integral (s^-1) |
|---|---:|
{table}

Separar as fases evita esconder uma reversão de sinal pela integral de toda a trajetória:

| Termo | Até o pico rastreado (s^-1) | Depois do pico (s^-1) |
|---|---:|---:|
{phase_table}

No ponto de entrada resolvido da parcela representativa, a vorticidade horizontal do perfil base interpolado é (xi,eta)=({environment_xi:.5g}, {environment_eta:.5g}) s^-1. A comparação com a vorticidade registrada informa quanto já estava presente ao entrar na região amostrada; não atribui retroativamente processos que ocorreram antes dessa entrada. As estatísticas das {result['parcel_count']} parcelas por fase estão em `parcel_phase_contributions.csv` e `summary.json`, para que a conclusão não dependa apenas da parcela representativa.

O rotacional da força de flutuabilidade produz as componentes horizontais, e elas foram preservadas em `lagrangian_budget.csv` como `buoyancy_xi` e `buoyancy_eta`. Seu termo vertical direto não deve ser confundido com a posterior inclinação dessas componentes. As contribuições de arrasto, LES, Coriolis, projeção e contornos também são calculadas a partir das mudanças efetivamente aplicadas. Não foi pressuposto que a vorticidade horizontal viesse toda do ambiente ou toda da baroclinicidade.

## O que permanece sugestivo ou não identificado

- A diferença contínuo–discreto da advecção não é desprezada: razão mediana de seu RMS para o RMS do rotacional advectivo={summary['continuum_mismatch_relative_median']:.3g}. O solver usa advecção conservativa em fluxo; em escoamento anelástico, a diferença para a forma material inclui termos de divergência, além de discretização, limitadores, recentramento e divisão de etapas. O campo `advection_operator_remainder` contém essa diferença; **não é automaticamente difusão numérica**.
- O orçamento material é uma integração interpolada, com residual explícito. Em {result['valid_budget_intervals']} de {result['total_budget_intervals']} intervalos-parcela válidos, RMS da mudança de zeta={delta:.3e} s^-1 e RMS do residual={rms:.3e} s^-1, razão={relative:.3g}. A precisão da atribuição física está limitada por esse residual, mesmo quando os incrementos Eulerianos fecham a arredondamento.
- Sensibilidade das posições: RK4 de 2 s versus 1 s, mediana {summary['trajectory_sensitivity']['rk4_2s_vs_1s_median_m']:.3g} m e máximo {summary['trajectory_sensitivity']['rk4_2s_vs_1s_max_m']:.3g} m; cadência nominal 30 s versus 60 s, mediana {summary['trajectory_sensitivity']['cadence_30_vs_60s_median_m']:.3g} m e máximo {summary['trajectory_sensitivity']['cadence_30_vs_60s_max_m']:.3g} m. As estatísticas consideram apenas pares finitos, cujas contagens estão no JSON; trajetórias censuradas não foram consideradas convergidas.
- A malha de 600 m e o primeiro centro acima do solo limitam a interpretação do vórtice e do funil ao nível resolvido. A inclinação pode transformar vorticidade horizontal preexistente em vertical, mas a origem anterior ao segmento amostrado não pode ser atribuída retroativamente.

## Encerramento

O produto separa a evolução de uma estrutura rotativa, a origem material resolvida, o orçamento discreto, o orçamento contínuo aproximado e a presença de condensado. Não aplica metas de intensidade, não insere vorticidade, não altera a termodinâmica e não propõe mudanças físicas. Onde a identidade, o orçamento ou a trajetória não sustentam uma atribuição, a limitação permanece explícita.
'''
    (out/'RELATORIO_FINAL.md').write_text(text,encoding='utf-8')


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument('sequence',type=Path)
    parser.add_argument('--out',type=Path); parser.add_argument('--seed-time',type=float,default=1800.)
    args=parser.parse_args(); out=args.out or args.sequence.parent/'analysis'; out.mkdir(parents=True,exist_ok=True)
    with h5py.File(args.sequence,'r') as f:
        if f.attrs['status']!='complete': raise RuntimeError('Analyze only a closed, complete sequence')
        summary,peak,seeds,track=audit_and_track(f,out,args.seed_time)
        print('RK4 trajectories',len(seeds),'seed index',peak,flush=True)
        positions=trajectories(f,peak,seeds,2.)
        fine=trajectories(f,peak,seeds,1.)
        coarse=trajectories(f,peak,seeds,2.,stride=2)
        distances=lambda other:np.concatenate([np.linalg.norm(positions[i]-p,axis=1) for i,p in other.items()])
        d=distances(fine); c=distances(coarse)
        summary['trajectory_sensitivity']=dict(rk4_2s_vs_1s_max_m=float(np.nanmax(d)),rk4_2s_vs_1s_median_m=float(np.nanmedian(d)),
                                               cadence_30_vs_60s_max_m=float(np.nanmax(c)),cadence_30_vs_60s_median_m=float(np.nanmedian(c)),
                                               finite_pairs_rk4=int(np.isfinite(d).sum()),finite_pairs_cadence=int(np.isfinite(c).sum()))
        result=lagrangian_budget(f,positions,out)
        summary['lagrangian']={k:v for k,v in result.items() if k not in ('records','budget')}
        summary['sequence_status']=f.attrs['status']; summary['metadata']=json.loads(f.attrs['metadata'])
        summary['analysis_source_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                          for p in (Path(__file__),ROOT/'src/storm_dynamics/diagnostic_capture.py')}
        mismatches=[]
        with zipfile.ZipFile(out/'source_snapshot.zip','w',compression=zipfile.ZIP_DEFLATED) as archive:
            for name,expected in summary['metadata'].get('source_sha256',{}).items():
                path=ROOT/name
                if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
                    mismatches.append(name); continue
                archive.write(path,name)
            if 'scripts/analyze_diagnostic_sequence.py' not in {n.replace('\\','/') for n in summary['metadata'].get('source_sha256',{})}:
                archive.write(Path(__file__),'scripts/analyze_diagnostic_sequence.py')
        summary['source_snapshot_mismatches']=mismatches
        write_report(summary,result,track,out,f)
        def json_safe(value):
            if isinstance(value,dict): return {k:json_safe(v) for k,v in value.items()}
            if isinstance(value,(list,tuple)): return [json_safe(v) for v in value]
            if isinstance(value,float) and not np.isfinite(value): return None
            return value
        (out/'summary.json').write_text(json.dumps(json_safe(summary),indent=2,default=str,allow_nan=False),encoding='utf-8')
        print(json.dumps({k:v for k,v in summary.items() if k!='metadata'},indent=2),flush=True)


if __name__=='__main__': main()
