"""Offline attribution of the STRONG-minus-WEAK vertical-vorticity contrast.

Uses only completed causal-pilot HDF5/CSV products.  The solver is never
imported or integrated.  Native velocity increments provide the exact
Eulerian operator closure; saved RK4 trajectories provide the approximate
material budget with its residual kept explicit.
"""
from __future__ import annotations

from pathlib import Path
import csv
import hashlib
import json
import math
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts")); sys.path.insert(0,str(ROOT/"src"))
from analyze_diagnostic_sequence import candidate_mask, native
from analyze_tilting_reversal import field_bundle
from storm_dynamics.diagnostic_capture import centered, curl, grad, kinematics

CASES={"WEAK-EVAP":ROOT/"outputs"/"causal_evap_weak",
       "CONTROL":ROOT/"outputs"/"causal_evap_control",
       "STRONG-EVAP":ROOT/"outputs"/"causal_evap_strong"}
COLORS={"WEAK-EVAP":"#2878b5","CONTROL":"#333333","STRONG-EVAP":"#d9534f"}
EVENTS={"2790":14,"2820":15,"2850":16,"2940":19,"3000":21}
WINDOWS={"2790-2850":(14,16),"2790-2940":(14,19),"2790-3000":(14,21)}
PHASES={"A_2790-2820":(14,15),"B_2820-2850":(15,16),
        "C_2850-2940":(16,19),"D_2940-3000":(19,21)}
TERMS=("stretching","tilting","dilatation","les","transport_advection_remainder",
       "projection_pressure","buoyancy_direct","drag","coriolis","other_operators","residual")


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig") as stream: return list(csv.DictReader(stream))


def write_csv(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w",newline="",encoding="utf-8") as stream:
        w=csv.DictWriter(stream,fieldnames=keys);w.writeheader();w.writerows(rows)


def f(row,key,default=np.nan):
    try:return float(row[key])
    except (KeyError,TypeError,ValueError):return default


def sample(array,coords,points):
    return RegularGridInterpolator(coords,array,bounds_error=False,fill_value=np.nan)(points)


def component_mask(zeta,z,track):
    labels,_=candidate_mask(zeta,z)
    pos=tuple(int(track[k]) for k in ("peak_i","peak_j","peak_k"))
    ident=int(labels[pos])
    if ident==0: raise RuntimeError(f"Tracked peak outside diagnosed component: {track}")
    return labels==ident


def weighted_mean(array,mask,volume):
    weights=volume*mask
    return float(np.sum(array*weights)/np.sum(weights))


def axis_tilt(mask,zeta,x,y,z):
    zs=[];xs=[];ys=[]
    xx,yy=np.meshgrid(x,y,indexing="ij")
    for k,height in enumerate(z):
        take=mask[:,:,k]
        if take.sum()<3:continue
        weight=np.maximum(zeta[:,:,k],0)*take
        if weight.sum()<=0:continue
        zs.append(height);xs.append(float((weight*xx).sum()/weight.sum()));ys.append(float((weight*yy).sum()/weight.sum()))
    if len(zs)<3:return np.nan
    sx=np.polyfit(zs,xs,1)[0];sy=np.polyfit(zs,ys,1)[0]
    return float(np.degrees(np.arctan(np.hypot(sx,sy))))


def canonical_budget(row):
    other=sum(f(row,key) for key in ("external_forcing_zeta","guard_zeta","initial_bcs_zeta",
        "predictor_bcs_zeta","projection_bcs_zeta","transport_microphysics_bcs_zeta"))
    return dict(observed_change=f(row,"delta_zeta"),stretching=f(row,"stretching_zeta"),
        tilting=f(row,"tilting_zeta"),dilatation=f(row,"vector_dilatation_zeta"),les=f(row,"les_zeta"),
        transport_advection_remainder=f(row,"advection_operator_remainder_zeta"),
        projection_pressure=f(row,"projection_zeta"),buoyancy_direct=f(row,"buoyancy_zeta"),
        drag=f(row,"surface_drag_zeta"),coriolis=f(row,"coriolis_zeta"),other_operators=other,
        residual=f(row,"residual_zeta"))


def load_tables(directory):
    track={int(r["index"]):r for r in read_csv(directory/"analysis"/"vortex_track.csv") if r.get("zeta_low_max_s")}
    geom={(int(r["index"]),int(r["parcel"])):r for r in read_csv(directory/"tilting_analysis"/"parcel_geometry.csv")}
    parcels={(int(r["index"]),int(r["parcel"])):r for r in read_csv(directory/"analysis"/"parcels.csv")}
    budget={(int(r["index"]),int(r["parcel"])):canonical_budget(r) for r in read_csv(directory/"analysis"/"lagrangian_budget.csv") if r["valid"].lower()=="true"}
    return dict(track=track,geometry=geom,parcels=parcels,budget=budget)


def diagnose_case(name,directory,out):
    tables=load_tables(directory); samples=[]; components=[]; euler=[]; maps={}
    with h5py.File(directory/"sequence.h5","r") as h5:
        groups=list(h5["snapshots"].values());nk=int(h5.attrs["budget_nz_with_halo"])
        x,y,z=(h5[f"grid/{n}"][:] for n in ("xc","yc","zc"));z=z[:nk]
        nphys=int(np.count_nonzero(z<=2000));z=z[:nphys]
        coords=(x,y,z);dx=float(np.diff(x).mean());dy=float(np.diff(y).mean())
        zf=h5["grid/zf"][:nphys+1];dz=np.diff(np.minimum(zf,2000.0));volume=np.broadcast_to((dx*dy*dz)[None,None,:],(len(x),len(y),nphys))
        rho=h5["base/rho0"][:nphys];rhow=h5["base/rho0_wface"][:nphys+1]
        xx,yy=np.meshgrid(x,y,indexing="ij");previous_zeta=None
        for i,group in enumerate(groups):
            vel=tuple(group[n][:,:,:nphys+(n=="w")] for n in ("u","v","w")); uc,vc,wc=centered(vel)
            kin=kinematics(vel,coords);xi,eta,zeta=kin["omega"]
            ux,uy,uz=grad(uc,coords);vx,vy,vz=grad(vc,coords);wx,wy,wz=grad(wc,coords)
            zx,zy,zz=grad(zeta,coords);transport_h=-(uc*zx+vc*zy);transport_v=-wc*zz
            div3=ux+vy+wz; ch=-(ux+vy)
            massdiv=(np.diff(vel[0],axis=0)*rho/float(np.diff(x).mean())+
                     np.diff(vel[1],axis=1)*rho/float(np.diff(y).mean())+
                     np.diff(vel[2]*rhow,axis=2)/np.diff(h5["grid/zf"][:nphys+1]))
            fields,_=field_bundle(group,h5,nk)
            fields={key:value[:,:,:nphys] for key,value in fields.items()}
            pd=group["p_dyn"][:,:,:nphys];pdx,pdy,pdz=grad(pd,coords)
            points=np.array([[f(tables["parcels"][(i,n)],k) for k in ("x_m","y_m","z_m")] for n in range(27)])
            arrays=dict(zeta=zeta,wz=wz,stretching=zeta*wz,convergence_h=ch,divergence_3d=div3,
                mass_divergence=massdiv,transport_horizontal=transport_h,transport_vertical=transport_v,
                omega_h=np.hypot(xi,eta),buoyancy_curl_h=np.hypot(fields["buoyancy_curl_xi"],fields["buoyancy_curl_eta"]),
                Bx=fields["Bx"],By=fields["By"],p_dyn=pd,p_dyn_gradient=np.sqrt(pdx*pdx+pdy*pdy+pdz*pdz))
            sampled={key:sample(value,coords,points) for key,value in arrays.items()}
            for n in range(27):samples.append(dict(case=name,index=i,parcel=n,time_s=float(group.attrs["time_s"]),
                **{key:float(value[n]) for key,value in sampled.items()}))
            track=tables["track"][i];mask=component_mask(zeta,z,track)
            k=int(track["peak_k"]);mask2=mask[:,:,k];area=mask2.sum()*dx*dy
            cx=f(track,"center_x_m");cy=f(track,"center_y_m");disk=np.hypot(xx-cx,yy-cy)<=3000
            gamma_component=float((zeta[:,:,k]*mask2).sum()*dx*dy)
            gamma_disk=float((zeta[:,:,k]*disk).sum()*dx*dy)
            wsearch=np.where(disk,wc[:,:,k],-np.inf);csearch=np.where(disk,ch[:,:,k],-np.inf)
            zw=np.unravel_index(np.argmax(wsearch),wsearch.shape);zc=np.unravel_index(np.argmax(csearch),csearch.shape)
            zp=(int(track["peak_i"]),int(track["peak_j"]))
            components.append(dict(case=name,index=i,time_s=float(group.attrs["time_s"]),
                zeta_max_s=f(track,"zeta_low_max_s"),zeta_volume_mean_s=weighted_mean(zeta,mask,volume),
                component_volume_km3=float(np.sum(volume*mask)/1e9),component_area_500m_km2=float(area/1e6),
                equivalent_radius_500m_km=float(np.sqrt(area/math.pi)/1000) if area else np.nan,
                circulation_component_500m_m2s=gamma_component,circulation_disk_3km_500m_m2s=gamma_disk,
                concentration_max_to_component_mean=float(f(track,"zeta_low_max_s")/(gamma_component/area)) if area and gamma_component else np.nan,
                axis_tilt_deg=axis_tilt(mask,zeta,x,y,z),z_min_m=f(track,"z_min_m"),z_max_m=f(track,"z_max_m"),
                zeta_to_wmax_distance_m=float(np.hypot(x[zp[0]]-x[zw[0]],y[zp[1]]-y[zw[1]])),
                zeta_to_convergence_max_distance_m=float(np.hypot(x[zp[0]]-x[zc[0]],y[zp[1]]-y[zc[1]])),
                p_dyn_core_minus_ring_Pa=f(track,"core_pressure_minus_ring_Pa"),
                p_dyn_gradient_disk_mean_Pam=float(np.mean(np.sqrt(pdx[:,:,k]**2+pdy[:,:,k]**2)[disk]))))
            if i in EVENTS.values():
                maps[i]=dict(time=float(group.attrs["time_s"]),x=x,y=y,z=z,zeta=zeta.copy(),w=wc.copy(),
                    convergence=ch.copy(),stretching=(zeta*wz).copy(),tilting=kin["tilting"].copy(),B=fields["B"].copy(),
                    buoyancy_curl_h=np.hypot(fields["buoyancy_curl_xi"],fields["buoyancy_curl_eta"]),p_dyn=pd.copy(),track=track)
            if i:
                ki=group["kinematic_integrals"]
                stretch=ki["stretching"][:,:,:nphys];tilt=ki["tilting"][:,:,:nphys]
                dil=ki["vector_dilatation"][2,:,:,:nphys];transport=ki["vector_advection"][2,:,:,:nphys]
                adv_native=curl(native(group["increments/advection"],nphys),coords)[2]
                remainder=adv_native-(stretch+tilt+dil+transport)
                stage={key:curl(native(value,nphys),coords)[2] for key,value in group["increments"].items() if key!="advection"}
                other=sum(stage.get(key,0) for key in ("external_forcing","guard","initial_bcs","predictor_bcs","projection_bcs","transport_microphysics_bcs"))
                budget=dict(stretching=stretch,tilting=tilt,dilatation=dil,transport_continuous=transport,
                    advection_operator_remainder=remainder,les=stage.get("les",np.zeros_like(zeta)),
                    projection_pressure=stage.get("projection",np.zeros_like(zeta)),
                    buoyancy_direct=stage.get("buoyancy",np.zeros_like(zeta)),drag=stage.get("surface_drag",np.zeros_like(zeta)),
                    coriolis=stage.get("coriolis",np.zeros_like(zeta)),other_operators=other)
                observed=zeta-previous_zeta;closure=observed-sum(budget.values())
                dt=float(group.attrs["time_s"])-float(groups[i-1].attrs["time_s"])
                peak_pos=(int(track["peak_i"]),int(track["peak_j"]),int(track["peak_k"]))
                euler.append(dict(case=name,index=i,time_s=float(group.attrs["time_s"]),dt_s=dt,
                    observed_change=weighted_mean(observed,mask,volume),closure=weighted_mean(closure,mask,volume),
                    **{key:weighted_mean(value,mask,volume) for key,value in budget.items()},
                    peak_observed_change=float(observed[peak_pos]),peak_closure=float(closure[peak_pos]),
                    **{"peak_"+key:float(value[peak_pos]) for key,value in budget.items()}))
            previous_zeta=zeta
    write_csv(out/f"parcel_samples_{name.lower().replace('-','_')}.csv",samples)
    write_csv(out/f"component_geometry_{name.lower().replace('-','_')}.csv",components)
    write_csv(out/f"eulerian_budget_{name.lower().replace('-','_')}.csv",euler)
    return dict(name=name,tables=tables,samples=samples,components=components,euler=euler,maps=maps)


def make_lagrangian_products(cases,out):
    interval=[];windows=[];phases=[];factor=[]
    bycase={name:{(r["index"],r["parcel"]):r for r in case["samples"]} for name,case in cases.items()}
    for name,case in cases.items():
        for (index,parcel),row in case["tables"]["budget"].items():
            if 15<=index<=21:
                dt=f(case["tables"]["geometry"][(index,parcel)],"time_s")-f(case["tables"]["geometry"][(index-1,parcel)],"time_s")
                interval.append(dict(case=name,index=index,parcel=parcel,time_s=f(case["tables"]["geometry"][(index,parcel)],"time_s"),dt_s=dt,
                    **row,**{key+"_rate_s2":value/dt for key,value in row.items()}))
    for label,(start,end) in {**WINDOWS,**PHASES}.items():
        kind="cumulative" if label in WINDOWS else "phase"
        for parcel in range(27):
            totals={}
            for name,case in cases.items():
                totals[name]={key:sum(case["tables"]["budget"][(i,parcel)][key] for i in range(start+1,end+1))
                              for key in ("observed_change",*TERMS)}
                z0=f(case["tables"]["geometry"][(start,parcel)],"zeta")
                ze=f(case["tables"]["geometry"][(end,parcel)],"zeta")
                totals[name].update(initial_zeta=z0,end_zeta=ze)
            for contrast,left,right in (("SW","STRONG-EVAP","WEAK-EVAP"),("SC","STRONG-EVAP","CONTROL"),("CW","CONTROL","WEAK-EVAP")):
                row=dict(window=label,kind=kind,parcel=parcel,contrast=contrast,start_index=start,end_index=end,
                    initial_zeta_contrast=totals[left]["initial_zeta"]-totals[right]["initial_zeta"],
                    endpoint_zeta_contrast=totals[left]["end_zeta"]-totals[right]["end_zeta"])
                row.update({key:totals[left][key]-totals[right][key] for key in ("observed_change",*TERMS)})
                (windows if kind=="cumulative" else phases).append(row)
    # Exact instantaneous product decomposition, using WEAK as the reference.
    for index in range(14,22):
        for parcel in range(27):
            w=bycase["WEAK-EVAP"][(index,parcel)];s=bycase["STRONG-EVAP"][(index,parcel)]
            dz=s["zeta"]-w["zeta"];dw=s["wz"]-w["wz"]
            product_delta=s["zeta"]*s["wz"]-w["zeta"]*w["wz"]
            factor.append(dict(index=index,parcel=parcel,time_s=s["time_s"],delta_stretching=product_delta,
                sampled_product_field_delta=s["stretching"]-w["stretching"],
                interpolation_product_mismatch=(s["stretching"]-w["stretching"])-product_delta,
                delta_zeta_times_wz_weak=dz*w["wz"],zeta_weak_times_delta_wz=w["zeta"]*dw,interaction=dz*dw,
                closure=product_delta-(dz*w["wz"]+w["zeta"]*dw+dz*dw)))
    write_csv(out/"lagrangian_interval_budgets.csv",interval);write_csv(out/"lagrangian_window_contrasts.csv",windows)
    write_csv(out/"lagrangian_phase_contrasts.csv",phases);write_csv(out/"stretching_factorization.csv",factor)
    return interval,windows,phases,factor


def summarize_ensemble(windows):
    rows=[]
    for window in WINDOWS:
        for contrast in ("SW","SC","CW"):
            selected=[r for r in windows if r["window"]==window and r["contrast"]==contrast]
            for key in ("initial_zeta_contrast","endpoint_zeta_contrast","observed_change",*TERMS):
                values=np.array([r[key] for r in selected])
                rows.append(dict(window=window,contrast=contrast,term=key,median=float(np.median(values)),minimum=float(values.min()),
                    maximum=float(values.max()),positive_count=int((values>0).sum()),negative_count=int((values<0).sum()),parcel_count=len(values)))
    return rows


def contrast_eulerian(cases,out):
    rows=[]
    lookup={name:{r["index"]:r for r in case["euler"]} for name,case in cases.items()}
    keys=[k for k in cases["CONTROL"]["euler"][0] if k not in ("case","index","time_s","dt_s")]
    for i in range(1,32):
        for label,left,right in (("SW","STRONG-EVAP","WEAK-EVAP"),("SC","STRONG-EVAP","CONTROL"),("CW","CONTROL","WEAK-EVAP")):
            a,b=lookup[left][i],lookup[right][i]
            rows.append(dict(index=i,time_s=a["time_s"],dt_s=a["dt_s"],contrast=label,
                **{key:a[key]-b[key] for key in keys}))
    aggregates=[]
    for window,(start,end) in WINDOWS.items():
        for contrast in ("SW","SC","CW"):
            chosen=[r for r in rows if r["contrast"]==contrast and start<r["index"]<=end]
            aggregates.append(dict(window=window,contrast=contrast,
                **{key:sum(r[key] for r in chosen) for key in keys}))
    write_csv(out/"eulerian_component_contrasts.csv",rows)
    write_csv(out/"eulerian_component_window_contrasts.csv",aggregates)
    return rows,aggregates


def figures(cases,interval,windows,ensemble,euler_contrasts,out):
    fd=out/"figures";fd.mkdir(exist_ok=True)
    def line(filename,source,key,ylabel,scale=1,central=False):
        fig,ax=plt.subplots(figsize=(11,5),constrained_layout=True)
        for name,case in cases.items():
            if source=="components":rows=case["components"]
            elif source=="samples":rows=[r for r in case["samples"] if r["parcel"]==13]
            else:rows=[r for r in interval if r["case"]==name and r["parcel"]==13]
            ax.plot([r["time_s"] for r in rows],[r[key]*scale for r in rows],label=name,color=COLORS[name])
        ax.axhline(0,color="gray",lw=.7);ax.set(xlabel="Tempo (s)",ylabel=ylabel,xlim=(2780,3010));ax.grid(alpha=.25);ax.legend();fig.savefig(fd/filename,dpi=160);plt.close(fig)
    line("01_zeta.png","components","zeta_max_s","ζ máxima do componente (s⁻¹)")
    line("02_Dzeta_Dt.png","interval","observed_change_rate_s2","Dζ/Dt material, parcela 13 (10⁻⁶ s⁻²)",1e6)
    line("03_stretching.png","interval","stretching_rate_s2","Stretching integrado / dt (10⁻⁶ s⁻²)",1e6)
    line("04_tilting.png","interval","tilting_rate_s2","Tilting integrado / dt (10⁻⁶ s⁻²)",1e6)
    line("05_les.png","interval","les_rate_s2","Curl LES integrado / dt (10⁻⁶ s⁻²)",1e6)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in case["samples"] if x["parcel"]==13 and 14<=x["index"]<=21]
        axs[0].plot([x["time_s"] for x in r],[x["transport_horizontal"]*1e6 for x in r],label=name,color=COLORS[name])
        axs[1].plot([x["time_s"] for x in r],[x["transport_vertical"]*1e6 for x in r],label=name,color=COLORS[name])
    axs[0].set_ylabel("−u∂xζ−v∂yζ (10⁻⁶ s⁻²)");axs[1].set_ylabel("−w∂zζ (10⁻⁶ s⁻²)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.axhline(0,color="gray",lw=.7);ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"06_transporte_horizontal_vertical.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in interval if x["case"]==name and x["parcel"]==13]
        axs[0].plot([x["time_s"] for x in r],[x["projection_pressure_rate_s2"]*1e12 for x in r],label=name,color=COLORS[name])
        p=[x for x in case["samples"] if x["parcel"]==13 and 14<=x["index"]<=21]
        axs[1].plot([x["time_s"] for x in p],[x["p_dyn_gradient"] for x in p],label=name,color=COLORS[name])
    axs[0].set_ylabel("Curl projeção (10⁻¹² s⁻²)");axs[1].set_ylabel("|∇p_dyn| (Pa m⁻¹)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"07_projection_pressao.png",dpi=160);plt.close(fig)
    sw=[r for r in interval if r["case"]=="STRONG-EVAP"]
    weak={(r["index"],r["parcel"]):r for r in interval if r["case"]=="WEAK-EVAP"}
    fig,ax=plt.subplots(figsize=(12,6),constrained_layout=True)
    for key in TERMS:
        vals=[];times=[]
        for r in sw:
            if r["parcel"]!=13:continue
            times.append(r["time_s"]);vals.append((r[key]-weak[(r["index"],13)][key])/r["dt_s"]*1e6)
        ax.plot(times,vals,label=key)
    ax.axhline(0,color="k",lw=.6);ax.set(xlabel="Tempo (s)",ylabel="Contraste SW de tendência (10⁻⁶ s⁻²)");ax.legend(ncol=2,fontsize=8);ax.grid(alpha=.25)
    fig.savefig(fd/"08_contrastes_todos_termos.png",dpi=160);plt.close(fig)
    central=[r for r in windows if r["contrast"]=="SW" and r["parcel"]==13]
    fig,ax=plt.subplots(figsize=(12,6),constrained_layout=True)
    for key in TERMS:ax.plot([int(r["window"].split("-")[1]) for r in central],[r[key]*1e6 for r in central],marker="o",label=key)
    ax.axhline(0,color="k",lw=.6);ax.set(xlabel="Fim da integral iniciada em 2790 s",ylabel="Contraste acumulado SW (10⁻⁶ s⁻¹)");ax.legend(ncol=2,fontsize=8);ax.grid(alpha=.25)
    fig.savefig(fd/"09_integrais_contrastes.png",dpi=160);plt.close(fig)
    target=next(r for r in central if r["window"]=="2790-3000")
    fig,axs=plt.subplots(1,2,figsize=(15,6),constrained_layout=True)
    for ax,use,title in ((axs[0],target,"Parcela central 13"),(axs[1],None,"Mediana das 27 parcelas")):
        values=[]
        for key in TERMS:
            if use is not None:values.append(use[key])
            else:values.append(next(r["median"] for r in ensemble if r["window"]=="2790-3000" and r["contrast"]=="SW" and r["term"]==key))
        ax.bar(np.arange(len(TERMS)),np.array(values)*1e6,color=["#4c78a8" if v>=0 else "#e45756" for v in values])
        ax.set_xticks(np.arange(len(TERMS)),TERMS,rotation=55,ha="right");ax.axhline(0,color="k",lw=.6);ax.set(ylabel="Contribuição SW (10⁻⁶ s⁻¹)",title=title)
    fig.suptitle("Waterfall contábil 2790–3000 s; residual mantido explícito");fig.savefig(fd/"10_waterfall_3000.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in case["samples"] if x["parcel"]==13 and 14<=x["index"]<=21]
        axs[0].plot([x["time_s"] for x in r],[x["wz"]*1e3 for x in r],label=name,color=COLORS[name])
        axs[1].plot([x["time_s"] for x in r],[x["convergence_h"]*1e3 for x in r],label=name,color=COLORS[name])
    axs[0].set_ylabel("∂w/∂z (10⁻³ s⁻¹)");axs[1].set_ylabel("Convergência horizontal (10⁻³ s⁻¹)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.axhline(0,color="gray",lw=.7);ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"11_wz_convergencia.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in case["components"] if 2780<=x["time_s"]<=3010]
        axs[0].plot([x["time_s"] for x in r],[x["circulation_component_500m_m2s"] for x in r],label=name,color=COLORS[name])
        axs[1].plot([x["time_s"] for x in r],[x["circulation_disk_3km_500m_m2s"] for x in r],label=name,color=COLORS[name])
    axs[0].set_ylabel("Γ no componente a ~500 m (m² s⁻¹)");axs[1].set_ylabel("Γ no disco de 3 km (m² s⁻¹)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"12_circulacao.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in case["components"] if 2780<=x["time_s"]<=3010]
        axs[0].plot([x["time_s"] for x in r],[x["component_area_500m_km2"] for x in r],label=name,color=COLORS[name])
        axs[1].plot([x["time_s"] for x in r],[x["component_volume_km3"] for x in r],label=name,color=COLORS[name])
    axs[0].set_ylabel("Área do componente a ~500 m (km²)");axs[1].set_ylabel("Volume do componente 0–2 km (km³)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"13_area_volume_componente.png",dpi=160);plt.close(fig)
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for name,case in cases.items():
        r=[x for x in case["samples"] if x["parcel"]==13 and 14<=x["index"]<=21]
        axs[0].plot([x["time_s"] for x in r],[x["omega_h"] for x in r],label=name,color=COLORS[name])
        axs[1].plot([x["time_s"] for x in r],[x["buoyancy_curl_h"]*1e6 for x in r],label=name,color=COLORS[name])
    axs[0].set_ylabel("|ωh| (s⁻¹)");axs[1].set_ylabel("|curl(B k)| (10⁻⁶ s⁻²)");axs[1].set_xlabel("Tempo (s)")
    for ax in axs:ax.grid(alpha=.25);ax.legend()
    fig.savefig(fd/"14_omega_h_gradientes_B.png",dpi=160);plt.close(fig)
    weak=cases["WEAK-EVAP"];strong=cases["STRONG-EVAP"]
    fig,axs=plt.subplots(2,5,figsize=(24,9),constrained_layout=True)
    for col,(label,index) in enumerate(EVENTS.items()):
        wm,sm=weak["maps"][index],strong["maps"][index];k=int(sm["track"]["peak_k"])
        for row,(key,unit,scale) in enumerate((("zeta","Δζ (10⁻⁴ s⁻¹)",1e4),("convergence","ΔCh (10⁻⁴ s⁻¹)",1e4))):
            delta=(sm[key][:,:,k]-wm[key][:,:,k])*scale;ax=axs[row,col]
            vmax=max(np.nanpercentile(abs(delta),99),1e-12);im=ax.pcolormesh(sm["x"]/1000,sm["y"]/1000,delta.T,cmap="RdBu_r",vmin=-vmax,vmax=vmax,shading="auto")
            cx=f(sm["track"],"center_x_m")/1000;cy=f(sm["track"],"center_y_m")/1000
            ax.scatter(cx,cy,c="gold",edgecolor="k",s=30);ax.set(xlim=(cx-8,cx+8),ylim=(cy-8,cy+8),title=f"{label} s: {unit}",xlabel="x (km)",ylabel="y (km)",aspect="equal");fig.colorbar(im,ax=ax)
    fig.savefig(fd/"15_mapas_horizontais_contraste.png",dpi=150);plt.close(fig)
    fig,axs=plt.subplots(2,5,figsize=(24,8),constrained_layout=True)
    for col,(label,index) in enumerate(EVENTS.items()):
        wm,sm=weak["maps"][index],strong["maps"][index];j=int(sm["track"]["peak_j"]);cx=f(sm["track"],"center_x_m")/1000
        for row,(key,unit,scale) in enumerate((("zeta","Δζ (10⁻⁴ s⁻¹)",1e4),("stretching","Δstretching (10⁻⁶ s⁻²)",1e6))):
            delta=(sm[key][:,j,:]-wm[key][:,j,:])*scale;ax=axs[row,col];vmax=max(np.nanpercentile(abs(delta),99),1e-12)
            im=ax.pcolormesh(sm["x"]/1000,sm["z"]/1000,delta.T,cmap="RdBu_r",vmin=-vmax,vmax=vmax,shading="auto")
            ax.set(xlim=(cx-8,cx+8),ylim=(0,2),title=f"{label} s: {unit}",xlabel="x (km)",ylabel="z (km)");fig.colorbar(im,ax=ax)
    fig.savefig(fd/"16_cortes_verticais_contraste.png",dpi=150);plt.close(fig)
    ens=[r for r in windows if r["window"]=="2790-3000" and r["contrast"]=="SW"]
    fig,ax=plt.subplots(figsize=(13,6),constrained_layout=True)
    vals=[[r[key]*1e6 for r in ens] for key in TERMS]
    ax.boxplot(vals,tick_labels=TERMS,showfliers=True);ax.axhline(0,color="k",lw=.6);ax.tick_params(axis="x",rotation=55);ax.set(ylabel="Contraste acumulado SW (10⁻⁶ s⁻¹)",title="27 parcelas vizinhas; não são réplicas independentes")
    fig.savefig(fd/"17_ensemble_27_parcelas.png",dpi=160);plt.close(fig)
    init=target["initial_zeta_contrast"]*1e6;end=target["endpoint_zeta_contrast"]*1e6
    fig,ax=plt.subplots(figsize=(14,5),constrained_layout=True);ax.axis("off")
    boxes=[("Vantagem já em 2790",f"Δζ inicial = {init:+.1f}×10⁻⁶ s⁻¹"),("2790–3000",f"tilt {target['tilting']*1e6:+.1f}; stretch {target['stretching']*1e6:+.1f}"),
           ("Termos restantes",f"LES {target['les']*1e6:+.1f}; dilat. {target['dilatation']*1e6:+.1f}; transp. rem. {target['transport_advection_remainder']*1e6:+.1f}"),
           ("Limite do orçamento",f"residual {target['residual']*1e6:+.1f}×10⁻⁶ s⁻¹"),("Em 3000",f"Δζ final = {end:+.1f}×10⁻⁶ s⁻¹")]
    for i,(title,text) in enumerate(boxes):
        x=.08+i*.21;ax.text(x,.5,title+"\n"+text,ha="center",va="center",bbox=dict(boxstyle="round,pad=.5",fc="#eef3f8",ec="#486581"))
        if i<4:ax.annotate("",xy=(x+.15,.5),xytext=(x+.075,.5),arrowprops=dict(arrowstyle="->",lw=1.6))
    ax.set_title("Mecanismo contábil: a diferença remanescente não excede o residual material")
    fig.savefig(fd/"18_mecanismo_compensatorio.png",dpi=170);plt.close(fig)


def main():
    out=ROOT/"outputs"/"zeta_compensation";out.mkdir(parents=True,exist_ok=True)
    cases={name:diagnose_case(name,path,out) for name,path in CASES.items()}
    interval,windows,phases,factor=make_lagrangian_products(cases,out)
    ensemble=summarize_ensemble(windows);write_csv(out/"ensemble_window_summary.csv",ensemble)
    euler_contrasts,euler_windows=contrast_eulerian(cases,out)
    figures(cases,interval,windows,ensemble,euler_contrasts,out)
    central=[r for r in windows if r["contrast"]=="SW" and r["parcel"]==13]
    counter=[]
    for r in central:
        for term in TERMS:
            counter.append(dict(window=r["window"],removed_term=term,actual_endpoint_contrast=r["endpoint_zeta_contrast"],
                diagnostic_endpoint_without_term=r["endpoint_zeta_contrast"]-r[term]))
    write_csv(out/"diagnostic_counterfactuals.csv",counter)
    summary=dict(input_only=True,no_new_simulation=True,critical_indices=EVENTS,windows=WINDOWS,phases=PHASES,
        central_SW={r["window"]:r for r in central},
        eulerian_component_SW={r["window"]:r for r in euler_windows if r["contrast"]=="SW"},
        ensemble_SW={window:{term:next(x for x in ensemble if x["window"]==window and x["contrast"]=="SW" and x["term"]==term)
                    for term in ("initial_zeta_contrast","endpoint_zeta_contrast","observed_change",*TERMS)} for window in WINDOWS},
        sequence_sha256={name:hashlib.sha256((path/"sequence.h5").read_bytes()).hexdigest() for name,path in CASES.items()},
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        limitations=dict(material_budget="30 s trajectory quadrature; explicit residual is attribution uncertainty",
            advection_remainder="native discrete advection curl minus independently reconstructed continuous transport/stretch/tilt/dilatation; not pure physical transport or numerical diffusion",
            circulation="area integral of reconstructed zeta on the ~492 m slice; disk radius 3 km provides a fixed-size comparison",
            dt="same archived CONTROL dt imposed; no adaptive-dt counterfactual exists",
            resolution="600 m horizontal grid and 30 s output cadence"))
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=True),encoding="utf-8")
    print(json.dumps(summary["central_SW"],indent=2),flush=True)


if __name__=="__main__":main()
