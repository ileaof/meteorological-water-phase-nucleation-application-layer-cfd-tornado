"""Cross-case analysis for the frozen rain-evaporation causal pilot.

Reads the three completed restart branches and their independently generated
component/parcel diagnostics.  It never calls the solver.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_tilting_reversal import field_bundle


CASES = {
    "WEAK-EVAP": ROOT / "outputs" / "causal_evap_weak",
    "CONTROL": ROOT / "outputs" / "causal_evap_control",
    "STRONG-EVAP": ROOT / "outputs" / "causal_evap_strong",
}
FACTORS = {"WEAK-EVAP": 0.95, "CONTROL": 1.0, "STRONG-EVAP": 1.05}
COLORS = {"WEAK-EVAP": "#2878b5", "CONTROL": "#333333", "STRONG-EVAP": "#d9534f"}


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def fnum(row, key, default=np.nan):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def truth(value):
    return str(value).lower() == "true"


def wrap_deg(value):
    return (value + 180.0) % 360.0 - 180.0


def nearest(rows, time):
    return min(rows, key=lambda row: abs(fnum(row, "time_s") - time))


def cumulative_processes(groups):
    names = sorted({name for group in groups[1:] for name in group["microphysics_integrals"].keys()})
    running = {name: 0.0 for name in names}
    rows = []
    for i, group in enumerate(groups):
        if i:
            for name in names:
                if name in group["microphysics_integrals"]:
                    running[name] += float(group["microphysics_integrals"][name][()])
        rows.append(dict(time_s=float(group.attrs["time_s"]), **running))
    return rows


def case_diagnostics(name, directory):
    track = read_csv(directory / "analysis" / "vortex_track.csv")
    track_valid = [row for row in track if row.get("zeta_low_max_s")]
    geometry_all = read_csv(directory / "tilting_analysis" / "parcel_geometry.csv")
    geometry = [row for row in geometry_all if int(row["parcel"]) == 13]
    lag_budget = [row for row in read_csv(directory / "analysis" / "lagrangian_budget.csv")
                  if int(row["parcel"]) == 13 and truth(row["valid"])]
    summary = json.loads((directory / "analysis" / "summary.json").read_text(encoding="utf-8"))
    tilt_summary = json.loads((directory / "tilting_analysis" / "summary.json").read_text(encoding="utf-8"))

    with h5py.File(directory / "sequence.h5", "r") as h5:
        groups = list(h5["snapshots"].values())
        nk = int(h5.attrs["budget_nz_with_halo"])
        x, y, z = (h5[f"grid/{axis}"][:] for axis in ("xc", "yc", "zc"))
        z = z[:nk]
        zf = h5["grid/zf"][:nk + 1]
        nphys = int(np.count_nonzero(z <= 2000.0))
        dx, dy = float(np.diff(x).mean()), float(np.diff(y).mean())
        ksurface = int(np.argmin(abs(z - 100.0)))
        kupdraft = int(np.argmin(abs(z - 500.0)))
        process_rows = cumulative_processes(groups)
        surface_rows = []
        precip = {}
        track_by_index = {int(row["index"]): row for row in track_valid}
        for i, group in enumerate(groups):
            fields, _ = field_bundle(group, h5, nk)
            dtv = fields["dtv_base"]
            B = fields["B"]
            tvs, Bs = dtv[:, :, ksurface], B[:, :, ksurface]
            cold = tvs < -1.0
            weights = np.where(cold, -tvs, 0.0)
            xx, yy = np.meshgrid(x, y, indexing="ij")
            cold_mass = weights.sum()
            cold_x = float((weights * xx).sum() / cold_mass) if cold_mass else np.nan
            cold_y = float((weights * yy).sum() / cold_mass) if cold_mass else np.nan
            contiguous = np.cumprod(dtv[:, :, :nphys] < -1.0, axis=2).astype(bool)
            clipped_faces = np.minimum(zf[:nphys + 1], 2000.0)
            depth = np.sum(contiguous * np.diff(clipped_faces)[None, None, :], axis=2)
            wup = fields["w"][:, :, kupdraft]
            up_weight = np.where(wup > 2.0, wup, 0.0)
            up_sum = up_weight.sum()
            up_x = float((up_weight * xx).sum() / up_sum) if up_sum else np.nan
            up_y = float((up_weight * yy).sum() / up_sum) if up_sum else np.nan
            tr = track_by_index.get(i, {})
            vortex_x, vortex_y = fnum(tr, "center_x_m"), fnum(tr, "center_y_m")
            gradB = np.hypot(fields["Bx"][:, :, ksurface], fields["By"][:, :, ksurface])
            row = dict(
                case=name, index=i, time_s=float(group.attrs["time_s"]),
                theta_v_area_lt_m05_km2=float(np.count_nonzero(tvs < -0.5) * dx * dy / 1e6),
                theta_v_area_lt_m1_km2=float(np.count_nonzero(cold) * dx * dy / 1e6),
                theta_v_area_lt_m2_km2=float(np.count_nonzero(tvs < -2.0) * dx * dy / 1e6),
                theta_v_min_K=float(tvs.min()),
                theta_v_min_0_2km_K=float(dtv[:, :, :nphys].min()),
                theta_v_cold_mean_K=float(tvs[cold].mean()) if cold.any() else np.nan,
                B_min_ms2=float(Bs.min()),
                B_min_0_2km_ms2=float(B[:, :, :nphys].min()),
                B_cold_mean_ms2=float(Bs[cold].mean()) if cold.any() else np.nan,
                cold_depth_max_m=float(depth.max()),
                cold_depth_mean_m=float(depth[cold].mean()) if cold.any() else np.nan,
                B_gradient_max_s2=float(gradB.max()),
                B_gradient_cold_mean_s2=float(gradB[cold].mean()) if cold.any() else np.nan,
                cold_centroid_x_m=cold_x, cold_centroid_y_m=cold_y,
                updraft_centroid_x_m=up_x, updraft_centroid_y_m=up_y,
                vortex_x_m=vortex_x, vortex_y_m=vortex_y,
                cold_to_updraft_distance_m=float(np.hypot(cold_x-up_x, cold_y-up_y)),
                vortex_to_updraft_distance_m=float(np.hypot(vortex_x-up_x, vortex_y-up_y)),
                vortex_to_cold_distance_m=float(np.hypot(vortex_x-cold_x, vortex_y-cold_y)),
                total_condensate_kgkg_sum=float(sum(group[q][:].sum() for q in ("ql","qi","qr","qs","qg","qh"))),
            )
            if i:
                for hyd, dataset in group["surface_precipitation_increment"].items():
                    precip[hyd] = precip.get(hyd, 0.0) + float(dataset[:].sum())
            row.update({f"surface_precip_{hyd}_kg": value for hyd, value in precip.items()})
            surface_rows.append(row)
        for j in range(1, len(surface_rows)-1):
            dr = (math.sqrt(surface_rows[j+1]["theta_v_area_lt_m1_km2"]*1e6/math.pi) -
                  math.sqrt(surface_rows[j-1]["theta_v_area_lt_m1_km2"]*1e6/math.pi))
            dt = surface_rows[j+1]["time_s"] - surface_rows[j-1]["time_s"]
            surface_rows[j]["equivalent_edge_speed_ms"] = dr/dt
        surface_rows[0]["equivalent_edge_speed_ms"] = np.nan
        surface_rows[-1]["equivalent_edge_speed_ms"] = np.nan
        metadata = json.loads(h5.attrs["metadata"])
        sequence_info = dict(
            sha256=hashlib.sha256((directory / "sequence.h5").read_bytes()).hexdigest(),
            bytes=(directory / "sequence.h5").stat().st_size,
            status=str(h5.attrs["status"]), frames=len(groups), state_nz=nk,
            start_s=float(groups[0].attrs["time_s"]), end_s=float(groups[-1].attrs["time_s"]),
            metadata=metadata,
        )

    first_low = next((row for row in track_valid if truth(row["touches_lowest_resolved_level"])), None)
    last_low = [row for row in track_valid if truth(row["touches_lowest_resolved_level"])][-1]
    material_peak = max(geometry, key=lambda row: fnum(row, "zeta"))
    crossing = next(row for row in tilt_summary["central_crossings"] if row["estimator"] == "T_product")
    weak_event = summary.get("sustained_20percent_weakening")
    cold_area_peak = max(surface_rows, key=lambda row: row["theta_v_area_lt_m1_km2"])
    cold_intensity_peak = min(surface_rows, key=lambda row: row["theta_v_min_K"])
    tz_min = min((row for row in geometry if fnum(row,"time_s") >= summary["peak"]["time_s"]),
                 key=lambda row: fnum(row,"T_product"))
    a, b = nearest(geometry, 2790.0), nearest(geometry, 2940.0)
    events = dict(
        A_intensification_onset_s=f"<={surface_rows[0]['time_s']:.6f} (censored by restart)",
        B_first_low_connection_s=fnum(first_low,"time_s") if first_low else None,
        C_low_zeta_peak_s=float(summary["peak"]["time_s"]),
        D_material_zeta_peak_s=fnum(material_peak,"time_s"),
        E_sustained_Tz_inversion_bracket_s=[float(crossing["t_positive_s"]),float(crossing["t_negative_s"])],
        F_persistent_zeta_loss_s=float(weak_event["time_s"]) if weak_event else None,
        G_last_low_connection_s=f">={fnum(last_low,'time_s'):.6f} (right-censored)",
        H_max_cold_area_s=float(cold_area_peak["time_s"]),
        H_max_cold_intensity_s=float(cold_intensity_peak["time_s"]),
    )
    metrics = dict(
        factor=FACTORS[name],
        rain_evaporation_cumulative_kg=process_rows[-1]["rain_evaporation_mass"],
        evaporative_cooling_cumulative_J=-process_rows[-1]["rain_evaporation_latent_energy"],
        cold_pool_area_max_km2=cold_area_peak["theta_v_area_lt_m1_km2"],
        theta_v_min_K=cold_intensity_peak["theta_v_min_K"],
        theta_v_min_0_2km_K=min(row["theta_v_min_0_2km_K"] for row in surface_rows),
        B_min_ms2=min(row["B_min_ms2"] for row in surface_rows),
        B_min_0_2km_ms2=min(row["B_min_0_2km_ms2"] for row in surface_rows),
        cold_depth_max_m=max(row["cold_depth_max_m"] for row in surface_rows),
        cold_depth_mean_at_low_zeta_peak_m=nearest(surface_rows,summary["peak"]["time_s"])["cold_depth_mean_m"],
        B_gradient_max_s2=max(row["B_gradient_max_s2"] for row in surface_rows),
        low_zeta_peak_time_s=float(summary["peak"]["time_s"]),
        low_zeta_peak_s=float(summary["peak"]["zeta_low_max_s"]),
        material_zeta_peak_time_s=fnum(material_peak,"time_s"),
        material_zeta_peak_s=fnum(material_peak,"zeta"),
        Tz_inversion_first_negative_s=float(crossing["t_negative_s"]),
        Tz_inversion_last_positive_s=float(crossing["t_positive_s"]),
        Tz_min_s2=fnum(tz_min,"T_product"),
        grad_w_rotation_2790_2940_deg=wrap_deg(fnum(b,"beta_deg")-fnum(a,"beta_deg")),
        omega_h_rotation_2790_2940_deg=wrap_deg(fnum(b,"alpha_deg")-fnum(a,"alpha_deg")),
        persistent_zeta_loss_s=float(weak_event["time_s"]) if weak_event else None,
        last_low_connection_s=fnum(last_low,"time_s"),
        low_connection_duration_lower_bound_s=fnum(last_low,"time_s")-fnum(first_low,"time_s") if first_low else None,
        low_connection_right_censored=True,
    )
    return dict(name=name, directory=directory, track=track_valid, geometry=geometry,
                geometry_all=geometry_all, lag_budget=lag_budget, process=process_rows,
                surface=surface_rows, summary=summary, tilt_summary=tilt_summary,
                events=events, metrics=metrics, sequence=sequence_info)


def state_contrasts(cases):
    handles = {name:h5py.File(case["directory"] / "sequence.h5", "r") for name,case in cases.items()}
    try:
        groups = {name:list(h5["snapshots"].values()) for name,h5 in handles.items()}
        nk = int(handles["CONTROL"].attrs["budget_nz_with_halo"])
        z = handles["CONTROL"]["grid/zc"][:nk]
        nphys = int(np.count_nonzero(z <= 2000.0))
        rows=[]
        for i in range(len(groups["CONTROL"])):
            bundles={name:field_bundle(gs[i],handles[name],nk)[0] for name,gs in groups.items()}
            raw={name:{field:gs[i][field][:,:,:nphys] for field in ("T","theta","qv","qr")}
                 for name,gs in groups.items()}
            raw={name:{**fields,"theta_v":bundles[name]["dtv_base"][:,:,:nphys],"B":bundles[name]["B"][:,:,:nphys]}
                 for name,fields in raw.items()}
            base=dict(index=i,time_s=float(groups["CONTROL"][i].attrs["time_s"]))
            for pair,left,right in (("strong_minus_weak","STRONG-EVAP","WEAK-EVAP"),
                                    ("weak_minus_control","WEAK-EVAP","CONTROL"),
                                    ("strong_minus_control","STRONG-EVAP","CONTROL")):
                for field in raw[left]:
                    delta=raw[left][field]-raw[right][field]
                    base[f"{pair}_{field}_mean"] = float(delta.mean())
                    base[f"{pair}_{field}_rms"] = float(np.sqrt(np.mean(delta*delta)))
                    base[f"{pair}_{field}_maxabs"] = float(np.abs(delta).max())
            rows.append(base)
        return rows
    finally:
        for h5 in handles.values(): h5.close()


def plot_lines(cases, filename, panels, xlabel="Tempo absoluto (s)"):
    fig, axes = plt.subplots(len(panels), 1, figsize=(11, 3.5*len(panels)), sharex=True,
                             constrained_layout=True)
    axes=np.atleast_1d(axes)
    for ax,(source,key,label,scale) in zip(axes,panels):
        for name,case in cases.items():
            rows=case[source]
            ax.plot([fnum(r,"time_s") for r in rows],[fnum(r,key)*scale for r in rows],
                    label=name,color=COLORS[name],lw=1.7)
        ax.set_ylabel(label); ax.grid(alpha=.25); ax.legend()
    axes[-1].set_xlabel(xlabel)
    fig.savefig(filename,dpi=160); plt.close(fig)


def plot_required(cases, contrasts, out):
    figures=out/"figures"; figures.mkdir(parents=True,exist_ok=True)
    plot_lines(cases,figures/"01_evaporacao_integrada.png",[("process","rain_evaporation_mass","Massa de chuva evaporada acumulada (10⁹ kg)",1e-9),
        ("process","rain_evaporation_latent_energy","Energia latente acumulada (10¹⁵ J)",1e-15)])
    plot_lines(cases,figures/"02_area_cold_pool.png",[("surface","theta_v_area_lt_m05_km2","Área θv′ < −0,5 K (km²)",1),
        ("surface","theta_v_area_lt_m1_km2","Área θv′ < −1 K (km²)",1),("surface","theta_v_area_lt_m2_km2","Área θv′ < −2 K (km²)",1)])
    plot_lines(cases,figures/"03_theta_v_minimo.png",[("surface","theta_v_min_K","Mínimo θv′ no nível de 121,7 m (K)",1),
        ("surface","theta_v_cold_mean_K","Média θv′ em θv′ < −1 K (K)",1)])
    plot_lines(cases,figures/"04_B_minimo.png",[("surface","B_min_ms2","B mínimo no nível de 121,7 m (m s⁻²)",1),
        ("surface","B_cold_mean_ms2","B médio em θv′ < −1 K (m s⁻²)",1)])
    plot_lines(cases,figures/"05_zeta_baixa.png",[("track","zeta_low_max_s","ζ máxima no componente abaixo de 500 m (s⁻¹)",1),
        ("geometry","zeta","ζ material, parcela 13 (s⁻¹)",1)])
    plot_lines(cases,figures/"06_Tz.png",[("geometry","T_product","Tz = ωh·∇hw, parcela 13 (10⁻⁶ s⁻²)",1e6),
        ("geometry","Tx_product","ξ ∂w/∂x (10⁻⁶ s⁻²)",1e6),("geometry","Ty_product","η ∂w/∂y (10⁻⁶ s⁻²)",1e6)])
    # Native interval-integrated stretching is divided by each actual output interval.
    fig,ax=plt.subplots(figsize=(11,4.5),constrained_layout=True)
    for name,case in cases.items():
        rr=case["lag_budget"]
        times=[fnum(r,"time_s") for r in rr]
        dt=np.diff([fnum(case["geometry"][0],"time_s"),*times])
        ax.plot(times,[fnum(r,"stretching_zeta")/d*1e6 for r,d in zip(rr,dt)],label=name,color=COLORS[name])
    ax.axhline(0,color="gray",lw=.7); ax.set(xlabel="Tempo absoluto (s)",ylabel="Stretching médio do intervalo (10⁻⁶ s⁻²)")
    ax.grid(alpha=.25); ax.legend(); fig.savefig(figures/"07_stretching.png",dpi=160); plt.close(fig)
    plot_lines(cases,figures/"08_angulo.png",[("geometry","angle_deg","Ângulo(ωh, ∇hw) (graus)",1)])
    plot_lines(cases,figures/"09_orientacao_omega_h.png",[("geometry","alpha_deg","Orientação de ωh (graus)",1)])
    plot_lines(cases,figures/"10_orientacao_grad_w.png",[("geometry","beta_deg","Orientação de ∇hw (graus)",1)])
    plot_lines(cases,figures/"11_componentes_grad_w.png",[("geometry","wx","∂w/∂x (10⁻³ s⁻¹)",1e3),
        ("geometry","wy","∂w/∂y (10⁻³ s⁻¹)",1e3)])

    # Horizontal maps at the low-level-zeta peak and first negative Tz sample.
    event_indices=[14,15]; event_names=["pico ζ baixo","primeiro Tz negativo"]
    fig,axes=plt.subplots(3,2,figsize=(13,16),constrained_layout=True)
    for row,(name,case) in enumerate(cases.items()):
        with h5py.File(case["directory"]/"sequence.h5","r") as h5:
            groups=list(h5["snapshots"].values()); nk=int(h5.attrs["budget_nz_with_halo"])
            x=h5["grid/xc"][:]/1000; y=h5["grid/yc"][:]/1000; z=h5["grid/zc"][:nk]
            for col,(idx,event_name) in enumerate(zip(event_indices,event_names)):
                fields,_=field_bundle(groups[idx],h5,nk); ks=int(np.argmin(abs(z-121.7))); kw=int(np.argmin(abs(z-500)))
                ax=axes[row,col]; im=ax.pcolormesh(x,y,fields["dtv_base"][:,:,ks].T,cmap="coolwarm",vmin=-4,vmax=4,shading="auto")
                ww=fields["w"][:,:,kw].T
                lev=[v for v in (-4,-2,0,2,5,10) if ww.min()<v<ww.max()]
                if lev: ax.contour(x,y,ww,levels=lev,colors="k",linewidths=.55)
                tr=next(r for r in case["track"] if int(r["index"])==idx)
                ax.scatter(fnum(tr,"center_x_m")/1000,fnum(tr,"center_y_m")/1000,c="gold",edgecolor="k",s=35)
                ax.set(title=f"{name}: {event_name}, t={float(groups[idx].attrs['time_s']):.1f} s",xlabel="x (km)",ylabel="y (km)",aspect="equal")
                fig.colorbar(im,ax=ax,label="θv′ a 121,7 m (K)")
    fig.savefig(figures/"12_mapas_horizontais_eventos.png",dpi=150); plt.close(fig)

    fig,axes=plt.subplots(1,3,figsize=(16,5.5),constrained_layout=True)
    for ax,(name,case) in zip(axes,cases.items()):
        idx=15; point=next(r for r in case["geometry"] if int(r["index"])==idx)
        with h5py.File(case["directory"]/"sequence.h5","r") as h5:
            groups=list(h5["snapshots"].values()); nk=int(h5.attrs["budget_nz_with_halo"])
            fields,_=field_bundle(groups[idx],h5,nk); x=h5["grid/xc"][:]/1000; y=h5["grid/yc"][:]; z=h5["grid/zc"][:nk]/1000
            j=int(np.argmin(abs(y-fnum(point,"y_m")))); im=ax.pcolormesh(x,z,fields["B"][:,j,:].T,cmap="coolwarm",vmin=-.15,vmax=.15,shading="auto")
            ww=fields["w"][:,j,:].T; lev=[v for v in (-4,-2,0,2,5,10) if ww.min()<v<ww.max()]
            if lev: ax.contour(x,z,ww,levels=lev,colors="k",linewidths=.7)
            ax.scatter(fnum(point,"x_m")/1000,fnum(point,"z_m")/1000,c="gold",edgecolor="k",s=35)
            ax.set(xlim=(fnum(point,"x_m")/1000-10,fnum(point,"x_m")/1000+10),ylim=(0,2),xlabel="x (km)",ylabel="z (km)",title=f"{name}, t={float(groups[idx].attrs['time_s']):.1f} s")
            fig.colorbar(im,ax=ax,label="B (m s⁻²)")
    fig.savefig(figures/"13_cortes_verticais_evento_inversao.png",dpi=160); plt.close(fig)

    fig=plt.figure(figsize=(16,5.5),constrained_layout=True)
    for n,(name,case) in enumerate(cases.items(),1):
        ax=fig.add_subplot(1,3,n,projection="3d")
        for parcel in range(27):
            rr=[r for r in case["geometry_all"] if int(r["parcel"])==parcel]
            ax.plot([fnum(r,"x_m")/1000 for r in rr],[fnum(r,"y_m")/1000 for r in rr],[fnum(r,"z_m")/1000 for r in rr],alpha=.25,color=COLORS[name])
        rr=case["geometry"]
        ax.scatter([fnum(r,"x_m")/1000 for r in rr],[fnum(r,"y_m")/1000 for r in rr],[fnum(r,"z_m")/1000 for r in rr],c=[fnum(r,"T_product") for r in rr],cmap="coolwarm",s=14)
        ax.set(title=name,xlabel="x (km)",ylabel="y (km)",zlabel="z (km)")
    fig.savefig(figures/"14_trajetorias_lagrangianas.png",dpi=160); plt.close(fig)

    fig,axes=plt.subplots(2,1,figsize=(11,8),constrained_layout=True)
    for name,case in cases.items():
        e=case["metrics"]["Tz_inversion_first_negative_s"]
        axes[0].plot([fnum(r,"time_s")-e for r in case["geometry"]],[fnum(r,"T_product")*1e6 for r in case["geometry"]],label=name,color=COLORS[name])
        p=case["metrics"]["low_zeta_peak_time_s"]
        axes[1].plot([fnum(r,"time_s")-p for r in case["track"]],[fnum(r,"zeta_low_max_s") for r in case["track"]],label=name,color=COLORS[name])
    axes[0].axvline(0,color="gray",ls=":"); axes[0].axhline(0,color="gray",lw=.6); axes[0].set(xlabel="Tempo relativo à inversão sustentada de Tz (s)",ylabel="Tz (10⁻⁶ s⁻²)")
    axes[1].axvline(0,color="gray",ls=":"); axes[1].set(xlabel="Tempo relativo ao pico de ζ baixa (s)",ylabel="ζ máxima abaixo de 500 m (s⁻¹)")
    for ax in axes: ax.grid(alpha=.25); ax.legend()
    fig.savefig(figures/"15_comparacao_alinhada_por_eventos.png",dpi=160); plt.close(fig)

    # Measured strong-minus-weak contrasts through the proposed mediation chain.
    weak,strong=cases["WEAK-EVAP"],cases["STRONG-EVAP"]
    at_peak=lambda case,key: nearest(case["surface"],case["metrics"]["low_zeta_peak_time_s"])[key]
    chain=[
        ("Evaporação\nacumulada",(strong["metrics"]["rain_evaporation_cumulative_kg"]-weak["metrics"]["rain_evaporation_cumulative_kg"])/1e9,"10⁹ kg"),
        ("Área fria\nno pico ζ",at_peak(strong,"theta_v_area_lt_m1_km2")-at_peak(weak,"theta_v_area_lt_m1_km2"),"km²"),
        ("|∇h w|\nparcela no pico",(fnum(nearest(strong["geometry"],2790),"grad_w_s")-fnum(nearest(weak["geometry"],2790),"grad_w_s"))*1e6,"10⁻⁶ s⁻¹"),
        ("Tz mínimo\npós-pico",(strong["metrics"]["Tz_min_s2"]-weak["metrics"]["Tz_min_s2"])*1e6,"10⁻⁶ s⁻²"),
        ("ζ baixa\nmáxima",(strong["metrics"]["low_zeta_peak_s"]-weak["metrics"]["low_zeta_peak_s"])*1e3,"10⁻³ s⁻¹"),
    ]
    fig,ax=plt.subplots(figsize=(14,4),constrained_layout=True); ax.axis("off")
    for i,(label,value,unit) in enumerate(chain):
        x=.08+i*.21
        ax.text(x,.52,f"{label}\nΔ(S−W)={value:+.4g} {unit}",ha="center",va="center",fontsize=10,
                bbox=dict(boxstyle="round,pad=.5",facecolor="#edf2f7",edgecolor="#486581"))
        if i<len(chain)-1: ax.annotate("",xy=(x+.145,.52),xytext=(x+.075,.52),arrowprops=dict(arrowstyle="->",lw=1.7))
    ax.set_title("Cadeia causal medida; contrastes não implicam ordenação no sentido de H1/H2")
    fig.savefig(figures/"16_cadeia_causal.png",dpi=170); plt.close(fig)

    plot_lines(cases,figures/"17_profundidade_gradiente_posicao_cold_pool.png",[("surface","cold_depth_max_m","Profundidade fria contígua máxima (m)",1),
        ("surface","B_gradient_max_s2","Máximo |∇hB| (10⁻⁶ s⁻²)",1e6),("surface","vortex_to_updraft_distance_m","Distância vórtice–updraft (km)",1e-3)])
    # Secondary microphysical effects at the end of the branch.
    names=[k for k in cases["CONTROL"]["process"][-1] if k.endswith("_mass") and k!="rain_evaporation_mass"]
    fig,ax=plt.subplots(figsize=(14,6),constrained_layout=True)
    y=np.arange(len(names)); height=.24
    ctrl=cases["CONTROL"]["process"][-1]
    for j,name in enumerate(("WEAK-EVAP","STRONG-EVAP")):
        vals=[]
        for key in names:
            denom=max(abs(ctrl[key]),1.0)
            vals.append((cases[name]["process"][-1][key]-ctrl[key])/denom*100)
        ax.barh(y+(j-.5)*height,vals,height=height,label=f"{name} − CONTROL",color=COLORS[name])
    ax.set_yticks(y,[k.removesuffix("_mass") for k in names]); ax.set_xlabel("Diferença acumulada relativa ao CONTROL (%)")
    ax.axvline(0,color="k",lw=.6); ax.legend(); fig.savefig(figures/"18_efeitos_microfisicos_secundarios.png",dpi=160); plt.close(fig)
    plot_lines(cases,figures/"19_velocidade_borda_equivalente.png",[("surface","equivalent_edge_speed_ms","Taxa de variação do raio de área equivalente (m s⁻¹)",1)])

    fig,axes=plt.subplots(3,2,figsize=(13,11),sharex=True,constrained_layout=True)
    fields=[("T","ΔT RMS (K)"),("theta","Δθ RMS (K)"),("qv","Δqv RMS (kg kg⁻¹)"),("qr","Δqr RMS (kg kg⁻¹)"),("theta_v","Δθv RMS (K)"),("B","ΔB RMS (m s⁻²)")]
    for ax,(field,label) in zip(axes.flat,fields):
        ax.plot([r["time_s"] for r in contrasts],[r[f"strong_minus_weak_{field}_rms"] for r in contrasts],color="#7b2cbf")
        ax.set_ylabel(label); ax.grid(alpha=.25)
    for ax in axes[-1]: ax.set_xlabel("Tempo absoluto (s)")
    fig.savefig(figures/"20_resposta_estado_strong_menos_weak.png",dpi=160); plt.close(fig)


def main():
    out=ROOT/"outputs"/"causal_evap_experiment"; out.mkdir(parents=True,exist_ok=True)
    cases={name:case_diagnostics(name,path) for name,path in CASES.items()}
    contrasts=state_contrasts(cases)
    for name,case in cases.items():
        write_csv(out/f"surface_{name.lower().replace('-','_')}.csv",case["surface"])
        write_csv(out/f"process_{name.lower().replace('-','_')}.csv",case["process"])
    write_csv(out/"state_contrasts.csv",contrasts)
    write_csv(out/"final_metrics.csv",[dict(case=name,**case["metrics"]) for name,case in cases.items()])
    secondary=[]
    for name,case in cases.items():
        for key,value in case["process"][-1].items():
            if key.endswith("_mass"):
                secondary.append(dict(case=name,process=key.removesuffix("_mass"),cumulative_mass_kg=value))
    write_csv(out/"secondary_microphysics.csv",secondary)
    plot_required(cases,contrasts,out)

    expected_frames=[]
    w,c,s=(cases[n] for n in ("WEAK-EVAP","CONTROL","STRONG-EVAP"))
    for i in range(len(c["surface"])):
        expected_frames.append(dict(
            time_s=c["surface"][i]["time_s"],
            evaporation_ordered=w["process"][i]["rain_evaporation_mass"] <= c["process"][i]["rain_evaporation_mass"] <= s["process"][i]["rain_evaporation_mass"],
            cold_area_ordered=w["surface"][i]["theta_v_area_lt_m1_km2"] <= c["surface"][i]["theta_v_area_lt_m1_km2"] <= s["surface"][i]["theta_v_area_lt_m1_km2"],
            cold_mean_strength_ordered=-w["surface"][i]["theta_v_cold_mean_K"] <= -c["surface"][i]["theta_v_cold_mean_K"] <= -s["surface"][i]["theta_v_cold_mean_K"],
            H1_Tz_ordered=fnum(s["geometry"][i],"T_product") <= fnum(c["geometry"][i],"T_product") <= fnum(w["geometry"][i],"T_product"),
            H1_zeta_ordered=fnum(s["track"][i],"zeta_low_max_s") <= fnum(c["track"][i],"zeta_low_max_s") <= fnum(w["track"][i],"zeta_low_max_s"),
        ))
    write_csv(out/"ordering_by_frame.csv",expected_frames)
    run_metadata={name:json.loads((case["directory"]/"metadata.json").read_text(encoding="utf-8"))
                  for name,case in cases.items()}
    stripped_configs=[]
    timestep_audit={}
    for name,case in cases.items():
        cfg=json.loads(json.dumps(run_metadata[name]["config"]))
        cfg["sim"]["physics"].pop("rain_evaporation_factor",None)
        stripped_configs.append(cfg)
        log=case["directory"]/"cfl_exceedances.log"; ratios=[]
        if log.exists():
            for line in log.read_text(encoding="utf-8").splitlines():
                _, imposed, branch_limit = map(float,line.split(","))
                ratios.append(imposed/branch_limit)
        timestep_audit[name]=dict(
            steps=run_metadata[name]["steps_executed"],
            archived_dt_min_s=run_metadata[name]["archived_dt_min_s"],
            archived_dt_max_s=run_metadata[name]["archived_dt_max_s"],
            imposed_dt_above_branch_limit_count=len(ratios),
            maximum_imposed_to_branch_limit_ratio=max(ratios,default=1.0),
            above_one_percent_count=sum(ratio>1.01 for ratio in ratios),
        )
    control_checks=run_metadata["CONTROL"]["control_validation"]
    result=dict(
        amplitude_selection=json.loads((ROOT/"outputs"/"causal_evap_screening_summary"/"selection.json").read_text(encoding="utf-8")),
        cases={name:{"events":case["events"],"metrics":case["metrics"],"sequence":{k:v for k,v in case["sequence"].items() if k!="metadata"},
                     "trajectory_sensitivity":case["summary"]["trajectory_sensitivity"],
                     "closure_max_s2":case["summary"]["closure_max_s2"],
                     "closure_relative_max":case["summary"]["closure_relative_max"]}
               for name,case in cases.items()},
        ordering_counts={key:int(sum(bool(row[key]) for row in expected_frames)) for key in expected_frames[0] if key!="time_s"},
        post_inversion_ordering_counts={
            key:int(sum(bool(row[key]) for row in expected_frames if 2820 <= row["time_s"] <= 3001))
            for key in ("H1_Tz_ordered","H1_zeta_ordered")
        },
        post_inversion_frame_count=int(sum(2820 <= row["time_s"] <= 3001 for row in expected_frames)),
        frames=len(expected_frames),
        state_contrast_maxima={key:max(float(row[key]) for row in contrasts) for key in contrasts[0] if key.endswith("_rms")},
        method=dict(
            intervention="restart at a common bitwise-identical archived state; factor active only from 2370.211552 s",
            synchronization="same archived CONTROL dt schedule and exact output times",
            cold_pool="theta_v perturbation from base state at z=121.66 m; thresholds -0.5/-1/-2 K",
            depth="contiguous from first stored cell upward where theta_v' < -1 K",
            edge_speed="time derivative of radius sqrt(area/pi), not a local gust-front speed",
            parcels="27 RK4 parcels around each equivalent tracked component's common-rule low-level peak; central ID 13 shown",
            tilting="independent velocity derivatives; Tz=xi*dw/dx+eta*dw/dy",
            event_resolution="30 s output brackets; identical event samples do not prove sub-output equality",
        ),
        execution_audit=dict(
            configurations_identical_after_removing_factor=stripped_configs[0]==stripped_configs[1]==stripped_configs[2],
            control_historical_bitwise_comparisons=len(control_checks),
            all_control_historical_comparisons_bitwise=all(row["bitwise"] and row["max_error"]==0 for row in control_checks),
            timestep=timestep_audit,
        ),
        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (out/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=True),encoding="utf-8")
    print(json.dumps(result,indent=2,allow_nan=True),flush=True)


if __name__=="__main__": main()
