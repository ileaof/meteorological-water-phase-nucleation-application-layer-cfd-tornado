"""Offline attribution of horizontal-vorticity provenance in the causal branches."""
from __future__ import annotations

import argparse, csv, hashlib, json
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASES = {
    "WEAK-EVAP": "weak",
    "CONTROL": "control",
    "STRONG-EVAP": "strong",
}
SOURCES = ("initial", "buoyancy", "surface_drag", "les", "coriolis",
           "projection", "boundary", "other", "advection_remainder")
COLORS = {s: c for s, c in zip(SOURCES, plt.cm.tab10.colors)}


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)


def sample_vector(a, coords, points):
    return np.column_stack([RegularGridInterpolator(coords, a[q], bounds_error=False,
        fill_value=np.nan)(points) for q in range(2)])


def case_analysis(root, case, stem):
    sim = root / f"vorticity_provenance_{stem}_staggered"
    ana = root / "vorticity_provenance_analysis" / stem
    parcel_rows = read_csv(ana / "tilting" / "parcel_geometry.csv")
    track_rows = {int(r["index"]): r for r in read_csv(ana / "base" / "vortex_track.csv")}
    by_index = {}
    for r in parcel_rows: by_index.setdefault(int(r["index"]), []).append(r)
    for v in by_index.values(): v.sort(key=lambda r: int(r["parcel"]))
    parcel_out, euler_out = [], []
    with h5py.File(sim / "provenance.h5") as f:
        x, y, z = (f[f"grid/{n}"][:] for n in ("xc", "yc", "zc"))
        coords = (x, y, z)
        keys = sorted(f["snapshots"])
        for idx, key in enumerate(keys):
            g = f["snapshots"][key]; rows = by_index[idx]
            pts = np.array([[float(r["x_m"]), float(r["y_m"]), float(r["z_m"])] for r in rows])
            grad = sample_vector(g["grad_w_h"][:], coords, pts)
            total = sample_vector(g["omega_total_h"][:], coords, pts)
            src = {s: sample_vector(g[f"omega_source_h/{s}"][:], coords, pts) for s in SOURCES}
            for n, r in enumerate(rows):
                rec = dict(case=case, index=idx, parcel=int(r["parcel"]), time_s=float(g.attrs["time_s"]),
                           x_m=pts[n,0], y_m=pts[n,1], z_m=pts[n,2],
                           xi=total[n,0], eta=total[n,1], wx=grad[n,0], wy=grad[n,1],
                           tilting=float(total[n] @ grad[n]), angle_deg=float(np.degrees(np.arccos(
                               np.clip((total[n] @ grad[n])/(np.linalg.norm(total[n])*np.linalg.norm(grad[n])+1e-300),-1,1)))))
                for s in SOURCES:
                    rec[f"xi_{s}"] = src[s][n,0]; rec[f"eta_{s}"] = src[s][n,1]
                    rec[f"tilting_{s}"] = float(src[s][n] @ grad[n])
                rec["tilting_reconstruction_error"] = rec["tilting"] - sum(rec[f"tilting_{s}"] for s in SOURCES)
                parcel_out.append(rec)
            # Eulerian disk centered on the tracked low-level component, 0--2 km.
            tr = track_rows[idx]; xx, yy = np.meshgrid(x, y, indexing="ij")
            mask2 = (xx-float(tr["center_x_m"]))**2 + (yy-float(tr["center_y_m"]))**2 <= 3000.0**2
            mask = mask2[:,:,None] & (z[None,None,:] <= 2000.)
            gt = g["grad_w_h"][:]; ot = g["omega_total_h"][:]; tt = np.sum(ot*gt, axis=0)
            er = dict(case=case,index=idx,time_s=float(g.attrs["time_s"]),cells=int(mask.sum()),
                      tilting_mean=float(np.mean(tt[mask])),tilting_positive_mean=float(np.mean(np.maximum(tt[mask],0))),
                      tilting_negative_mean=float(np.mean(np.minimum(tt[mask],0))),
                      omega_condition_index=float(g.attrs["omega_condition_index"]),
                      tilting_condition_index=float(g.attrs["tilting_condition_index"]),
                      tilting_advection_remainder_fraction=float(g.attrs["tilting_advection_remainder_fraction"]))
            for s in SOURCES:
                ts = np.sum(g[f"omega_source_h/{s}"][:]*gt,axis=0)
                er[f"tilting_{s}_mean"] = float(np.mean(ts[mask]))
            er["tilting_reconstruction_error"] = er["tilting_mean"]-sum(er[f"tilting_{s}_mean"] for s in SOURCES)
            euler_out.append(er)
    return parcel_out, euler_out


def summarize(parcels, euler):
    phases={"before":(2400,2760),"peak":(2760,2821),"after":(2821,3001),"late":(3001,3301)}
    out=[]
    for case in CASES:
        for phase,(a,b) in phases.items():
            p=[r for r in parcels if r["case"]==case and a<=r["time_s"]<b]
            e=[r for r in euler if r["case"]==case and a<=r["time_s"]<b]
            rec=dict(case=case,phase=phase,n_parcel_samples=len(p),n_euler_frames=len(e),
                     tilting_parcel_mean=float(np.nanmean([r["tilting"] for r in p])),
                     angle_parcel_mean_deg=float(np.nanmean([r["angle_deg"] for r in p])),
                     tilting_euler_mean=float(np.mean([r["tilting_mean"] for r in e])))
            for s in SOURCES:
                rec[f"tilting_{s}_parcel_mean"]=float(np.nanmean([r[f"tilting_{s}"] for r in p]))
                rec[f"tilting_{s}_euler_mean"]=float(np.mean([r[f"tilting_{s}_mean"] for r in e]))
            out.append(rec)
    return out


def contrasts(phases):
    rows=[]
    for phase in ("before","peak","after","late"):
        d={(r["case"],r["phase"]):r for r in phases}
        for label,a,b in (("SW","STRONG-EVAP","WEAK-EVAP"),("SC","STRONG-EVAP","CONTROL"),("CW","CONTROL","WEAK-EVAP")):
            rec=dict(contrast=label,phase=phase)
            for key in ("tilting_parcel_mean","angle_parcel_mean_deg","tilting_euler_mean")+tuple(f"tilting_{s}_parcel_mean" for s in SOURCES)+tuple(f"tilting_{s}_euler_mean" for s in SOURCES):
                rec[key]=float(d[(a,phase)][key])-float(d[(b,phase)][key])
            rows.append(rec)
    return rows


def figures(out, parcels, euler, phases):
    fd=out/"figures"; fd.mkdir(parents=True,exist_ok=True)
    # 27-parcel mean source contributions.
    fig,axs=plt.subplots(3,1,figsize=(12,12),sharex=True,constrained_layout=True)
    for ax,case in zip(axs,CASES):
        times=sorted(set(r["time_s"] for r in parcels if r["case"]==case))
        for s in SOURCES:
            vals=[np.nanmean([r[f"tilting_{s}"] for r in parcels if r["case"]==case and r["time_s"]==t]) for t in times]
            ax.plot(times,np.array(vals)*1e6,label=s,color=COLORS[s])
        ax.axhline(0,color="k",lw=.6); ax.axvline(2790,color="gray",ls=":"); ax.set_ylabel(case+"\n$T_z$ (10$^{-6}$ s$^{-2}$)")
    axs[-1].set_xlabel("Tempo (s)"); axs[0].legend(ncol=3,fontsize=8)
    fig.savefig(fd/"01_parcelas_fontes_tilting.png",dpi=170); plt.close(fig)
    # Total geometry and causal contrast.
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True,constrained_layout=True)
    for case in CASES:
        rs=[r for r in parcels if r["case"]==case and r["parcel"]==13]
        axs[0].plot([r["time_s"] for r in rs],[r["tilting"]*1e6 for r in rs],label=case)
        axs[1].plot([r["time_s"] for r in rs],[r["angle_deg"] for r in rs],label=case)
    axs[0].axhline(0,color="k",lw=.6); axs[1].axhline(90,color="k",lw=.6); axs[0].set_ylabel("Tilting parcela 13 (10$^{-6}$ s$^{-2}$)")
    axs[1].set(ylabel="Ângulo (graus)",xlabel="Tempo (s)"); axs[0].legend(); fig.savefig(fd/"02_geometria_parcela13.png",dpi=170); plt.close(fig)
    # Eulerian source contributions.
    fig,axs=plt.subplots(3,1,figsize=(12,12),sharex=True,constrained_layout=True)
    for ax,case in zip(axs,CASES):
        rs=[r for r in euler if r["case"]==case]
        for s in SOURCES: ax.plot([r["time_s"] for r in rs],[r[f"tilting_{s}_mean"]*1e6 for r in rs],label=s,color=COLORS[s])
        ax.axhline(0,color="k",lw=.6); ax.set_ylabel(case+"\nmean $T_z$ (10$^{-6}$ s$^{-2}$)")
    axs[-1].set_xlabel("Tempo (s)"); axs[0].legend(ncol=3,fontsize=8); fig.savefig(fd/"03_euleriano_fontes_tilting.png",dpi=170); plt.close(fig)
    # Phase table heatmap.
    arr=np.array([[next(r for r in phases if r["case"]==c and r["phase"]==p)[f"tilting_{s}_parcel_mean"] for s in SOURCES]
                  for c in CASES for p in ("before","peak","after","late")])*1e6
    fig,ax=plt.subplots(figsize=(12,8),constrained_layout=True); lim=max(abs(arr.min()),abs(arr.max()))
    im=ax.imshow(arr,aspect="auto",cmap="RdBu_r",vmin=-lim,vmax=lim)
    ax.set_xticks(range(len(SOURCES)),SOURCES,rotation=45,ha="right"); ax.set_yticks(range(len(arr)),[f"{c} {p}" for c in CASES for p in ("before","peak","after","late")])
    fig.colorbar(im,ax=ax,label="Tilting médio das 27 parcelas (10$^{-6}$ s$^{-2}$)"); fig.savefig(fd/"04_fases_fontes.png",dpi=170); plt.close(fig)
    # 3-D view of all fixed parcel IDs through the reversal window.
    fig=plt.figure(figsize=(11,8)); ax=fig.add_subplot(111,projection="3d")
    rs=[r for r in parcels if r["case"]=="CONTROL" and 2700<=r["time_s"]<=3001]
    sc=ax.scatter(np.array([r["x_m"] for r in rs])/1000,np.array([r["y_m"] for r in rs])/1000,
                  np.array([r["z_m"] for r in rs])/1000,c=np.array([r["tilting"] for r in rs])*1e6,
                  cmap="RdBu_r",s=12,vmin=-80,vmax=80)
    ax.set(xlabel="x (km)",ylabel="y (km)",zlabel="z (km)",title="CONTROL: 27 trajetórias, 2700–3000 s")
    fig.colorbar(sc,ax=ax,label="$T_z$ (10$^{-6}$ s$^{-2}$)"); fig.savefig(fd/"05_trajetorias_3d.png",dpi=170); plt.close(fig)


def spatial_maps(outputs,out):
    fd=out/"figures"; sim=outputs/"vorticity_provenance_control_staggered"/"provenance.h5"
    with h5py.File(sim) as f:
        x=f["grid/xc"][:]/1000; y=f["grid/yc"][:]/1000; z=f["grid/zc"][:]; k=int(np.argmin(abs(z-500)))
        groups=list(f["snapshots"].values()); targets=(2790,2820,2940,3000)
        chosen=[min(groups,key=lambda g:abs(float(g.attrs["time_s"])-t)) for t in targets]
        names=("total","buoyancy","projection","surface_drag","les")
        vals=[]
        for g in chosen:
            grad=g["grad_w_h"][:,:,:,k]
            vals.append([np.sum((g["omega_total_h"][:] if n=="total" else g[f"omega_source_h/{n}"][:])[:,:,:,k]*grad,axis=0)*1e6 for n in names])
        lim=np.nanpercentile(np.abs(np.asarray(vals)),99)
        fig,axs=plt.subplots(len(chosen),len(names),figsize=(18,14),constrained_layout=True)
        for i,g in enumerate(chosen):
            for j,n in enumerate(names):
                im=axs[i,j].pcolormesh(x,y,vals[i][j].T,cmap="RdBu_r",vmin=-lim,vmax=lim,shading="auto")
                axs[i,j].set_aspect("equal"); axs[i,j].set_title(f"{n}; {float(g.attrs['time_s']):.1f} s")
                if i==len(chosen)-1: axs[i,j].set_xlabel("x (km)")
                if j==0: axs[i,j].set_ylabel("y (km)")
        fig.colorbar(im,ax=axs,label="$T_z$ a ~500 m (10$^{-6}$ s$^{-2}$)",shrink=.7)
        fig.savefig(fd/"06_mapas_proveniencia_500m.png",dpi=160); plt.close(fig)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--outputs",type=Path,default=Path("outputs")); ap.add_argument("--out",type=Path,default=Path("outputs/vorticity_provenance_analysis")); a=ap.parse_args()
    a.out.mkdir(parents=True,exist_ok=True); parcels=[]; euler=[]
    for case,stem in CASES.items():
        p,e=case_analysis(a.outputs,case,stem); parcels+=p; euler+=e; print(case,len(p),len(e),flush=True)
    phases=summarize(parcels,euler); delta=contrasts(phases); write_csv(a.out/"parcel_source_attribution.csv",parcels); write_csv(a.out/"eulerian_source_attribution.csv",euler); write_csv(a.out/"phase_source_summary.csv",phases); write_csv(a.out/"causal_contrasts.csv",delta)
    figures(a.out,parcels,euler,phases)
    spatial_maps(a.outputs,a.out)
    closure={"max_parcel_abs":max(abs(r["tilting_reconstruction_error"]) for r in parcels),"max_euler_abs":max(abs(r["tilting_reconstruction_error"]) for r in euler),"max_saved_remainder_fraction":max(r["tilting_advection_remainder_fraction"] for r in euler)}
    meta={"cases":list(CASES),"sources":list(SOURCES),"closure":closure,"method":"trilinear sampling at unchanged 27 parcel positions; Eulerian 3-km disk about unchanged overlap-tracked center, z<=2 km","script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (a.out/"provenance_summary.json").write_text(json.dumps(meta,indent=2),encoding="utf-8")
    print(json.dumps(meta,indent=2),flush=True)

if __name__=="__main__": main()
