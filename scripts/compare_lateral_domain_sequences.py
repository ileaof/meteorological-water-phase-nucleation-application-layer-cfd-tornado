"""Compare 72 and 120 km sequences using common physical diagnostics."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_tornadogenesis_mechanism import (
    centered, component_mask, distance2, kinematics, periodic_delta,
    radial_metrics, track_axis,
)

RADIUS_M = 4200.0
TOP_M = 2000.0
PRIMARY = {
    "zeta_max_s-1": 1,
    "circulation_4200_m2_s": 1,
    "vtheta_max_m_s": 1,
    "zeta_halfmax_width_m": -1,
    "axis_tilt_surface_to_2km_m": -1,
}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def read_sequence(path):
    rows=[]; previous_surface=None
    with h5py.File(path,"r") as f:
        x,y,zfull,zf=(f[f"grid/{n}"][:] for n in ("xc","yc","zc","zf"))
        nk=int(np.searchsorted(zfull,TOP_M,side="right")); z=zfull[:nk]
        dx_m=float(x[1]-x[0]);dy_m=float(y[1]-y[0])
        base_tv=f["base/theta0"][:]*(1+.61*f["base/qv0"][:])
        for key,g in f["snapshots"].items():
            t=float(g.attrs["time_s"])
            if not 2790<=t<=3301: continue
            u,v,w=centered(g,nk);kin=kinematics(u,v,w,x,y,z);zeta=kin["zeta"]
            axis,previous_surface=track_axis(zeta,x,y,z,previous_surface)
            k=1;center=tuple(axis[k]);local=distance2(x,y,center)<=RADIUS_M**2
            gamma,_,_,rmw,vmax,pdef,radius,_=radial_metrics(
                zeta,u,v,np.asarray(g["p_dyn"][:,:,:nk]),x,y,center,k,dx_m,600.
            )
            z2=zeta[:,:,k];conv2=kin["convergence"][:,:,k]
            p2=np.asarray(g["p_dyn"][:,:,k]);ann=(radius>=4800)&(radius<=7200)
            panom=p2-float(p2[ann].mean())
            def diameter(mask):
                return float(2*np.sqrt(mask.sum()*dx_m*dy_m/np.pi))
            zhalf=local&(z2>=.5*z2[local].max())
            chalf=local&(conv2>=.5*conv2[local].max())
            phalf=local&(panom<=.5*pdef)
            high=component_mask(zeta,x,y,z,axis,.005)
            disp=np.hypot(periodic_delta(axis[:,0],axis[0,0],x[-1]-x[0]+dx_m),
                          periodic_delta(axis[:,1],axis[0,1],y[-1]-y[0]+dy_m))
            dz2=np.maximum(0,np.minimum(zf[1:nk+1],TOP_M)-np.maximum(zf[:nk],0))[None,None,:]
            volume=dx_m*dy_m*dz2
            cylinder=np.stack([distance2(x,y,axis[kk])<=RADIUS_M**2 for kk in range(nk)],axis=2)
            condensate=sum(np.asarray(g[n][:,:,:nk]) for n in ("ql","qi","qr","qs","qg","qh"))
            tv=np.asarray(g["theta"][:,:,k])*(1+.61*np.asarray(g["qv"][:,:,k])-condensate[:,:,k])
            cold=tv-base_tv[k] < -1
            rows.append({
                "time_s":t,"center_x_m":center[0],"center_y_m":center[1],
                "zeta_max_s-1":float(z2[local].max()),
                "zeta_max_0_2km_s-1":float(zeta[cylinder].max()),
                "zeta_signed_0_2km_m3_s":float((zeta*volume)[cylinder].sum()),
                "zeta_absolute_0_2km_m3_s":float((abs(zeta)*volume)[cylinder].sum()),
                "circulation_4200_m2_s":gamma[4],"vtheta_max_m_s":vmax,"rmw_m":rmw,
                "zeta_halfmax_width_m":diameter(zhalf),
                "convergence_halfmax_width_m":diameter(chalf),
                "pressure_halfdeficit_width_m":diameter(phalf),
                "pressure_dynamic_deficit_Pa":pdef,
                "convergence_max_s-1":float(conv2[local].max()),
                "stretching_conditional_mean_s-2":float(kin["stretching"][high].mean()) if high.any() else float("nan"),
                "axis_tilt_surface_to_2km_m":float(disp[-1]),"axis_tilt_max_m":float(disp.max()),
                "cold_pool_area_km2":float(cold.sum()*dx_m*dy_m/1e6),
                "cold_pool_fraction_near_vortex":float((cold&local).sum()/max(local.sum(),1)),
                "wmax_m_s":float(w[:,:,k][local].max()),
            })
        steps=np.asarray(f["steps"][:])
        geometry={"nx":len(x),"ny":len(y),"nz":len(zfull),"Lx_m":float(f["grid/xf"][-1]-f["grid/xf"][0]),
                  "Ly_m":float(f["grid/yf"][-1]-f["grid/yf"][0]),"Lz_m":float(zf[-1]-zf[0]),
                  "dx_m":dx_m,"dy_m":dy_m,"analysis_top_m":float(z[-1])}
    return rows,geometry,steps


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--control",type=Path,default=ROOT/"outputs/diagnostic_sequence_20260905/sequence.h5")
    ap.add_argument("--extended",type=Path,default=ROOT/"outputs/domain_extent_120km_600m_20260910/sequence.h5")
    ap.add_argument("--out",type=Path,default=ROOT/"outputs/domain_extent_120km_comparison_20260910")
    ap.add_argument("--self-test-identical",action="store_true")
    args=ap.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    control_rows,control_geometry,control_steps=read_sequence(args.control)
    extended_rows,extended_geometry,extended_steps=read_sequence(args.extended)
    data={72:control_rows,120:extended_rows}
    assert control_geometry["Lx_m"]==control_geometry["Ly_m"]==72000.,control_geometry
    expected_extended=72000. if args.self_test_identical else 120000.
    assert extended_geometry["Lx_m"]==extended_geometry["Ly_m"]==expected_extended,extended_geometry
    for geometry in (control_geometry,extended_geometry):
        assert geometry["dx_m"]==geometry["dy_m"]==600.,geometry
        assert geometry["Lz_m"]==15000. and geometry["nz"]==48,geometry
    times={d:np.array([r["time_s"] for r in rows]) for d,rows in data.items()}
    first=max(t[0] for t in times.values());last=min(t[-1] for t in times.values())
    common=np.unique(np.round(np.r_[np.arange(np.ceil(first/30)*30,last-1e-9,30),last],9))
    damping={}
    for domain,steps in ((72,control_steps),(120,extended_steps)):
        overlap=np.maximum(0,np.minimum(steps[:,0]+steps[:,1],last)-np.maximum(steps[:,0],first))
        active=overlap>0;calls=float(np.sum(4*overlap[active]/steps[active,1]))
        damping[str(domain)]={"overlapping_native_steps":int(active.sum()),
                              "fractional_operator_calls":calls,
                              "nominal_lowest_face_amplitude_decay_rate_s-1":float(-np.log(.95)*calls/(last-first))}
    excluded={"time_s","center_x_m","center_y_m"}
    metrics=sorted((set(data[72][0])&set(data[120][0]))-excluded)
    paired=[];stats=[]
    for metric in metrics:
        vals={d:np.interp(common,times[d],[float(r[metric]) for r in data[d]]) for d in data}
        a,b=vals[72],vals[120];valid=np.isfinite(a)&np.isfinite(b);delta=b-a
        scale=np.nanmedian(np.abs(a[valid]));fractional=float(np.nanmedian(delta[valid])/scale) if scale else float("nan")
        row={"metric":metric,"median_72km":float(np.nanmedian(a)),"median_120km":float(np.nanmedian(b)),
             "median_delta":float(np.nanmedian(delta)),"median_fractional_delta":fractional,
             "fraction_120km_greater":float(np.mean(b[valid]>a[valid])),
             "fraction_120km_less":float(np.mean(b[valid]<a[valid])),"valid_pairs":int(valid.sum())}
        stats.append(row)
        paired.extend({"time_s":float(t),"metric":metric,"value_72km":float(x),
                       "value_120km":float(y),"delta":float(y-x)} for t,x,y in zip(common,a,b))
    native=[]
    for domain,rows in data.items():
        native.extend({"domain_km":domain,**row} for row in rows)
    write_csv(args.out/"native_diagnostics.csv",native)
    write_csv(args.out/"paired_timeseries.csv",paired);write_csv(args.out/"metric_comparison.csv",stats)
    primary=[]
    for metric,direction in PRIMARY.items():
        r=next(x for x in stats if x["metric"]==metric)
        favorable=(r["fraction_120km_greater"] if direction>0 else r["fraction_120km_less"])
        primary.append({**r,"favorable_direction":direction,"favorable_pairs":int(round(favorable*r["valid_pairs"])),
                        "coherent_12_of_18":bool(favorable*r["valid_pairs"]>=12),
                        "material_20_percent":bool(direction*r["median_fractional_delta"]>=.20)})
    materially_favorable=sum(r["coherent_12_of_18"] and r["material_20_percent"] for r in primary)
    result={"status":"complete","common_times_s":common.tolist(),
            "geometry":{"control":control_geometry,"extended":extended_geometry},
            "nominal_top_damping":damping,
            "primary_metrics":primary,"all_metrics":stats,
            "screening_classification":("material favorable domain effect" if materially_favorable>=3 else
                                         "no dominant favorable domain effect" if materially_favorable<=1 else "mixed domain effect"),
            "materially_favorable_primary_metrics":materially_favorable,
            "limitations":["Single deterministic pair; not an ensemble uncertainty estimate.",
                           "Scalar diagnostics are interpolated; fields are not interpolated.",
                           "This comparison isolates lateral extent, not top height or damping law.",
                           "Adaptive trajectories may produce different step counts; nominal damping exposure is reported."]}
    (args.out/"summary.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    artifacts={}
    for path in sorted(args.out.iterdir()):
        if path.is_file() and path.name!="manifest.json":
            artifacts[path.name]={"bytes":path.stat().st_size,
                                  "sha256":hashlib.sha256(path.read_bytes()).hexdigest()}
    (args.out/"manifest.json").write_text(json.dumps({"artifacts":artifacts},indent=2),encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__": main()
