"""Offline, non-causal audit of the timestep-dependent top damping operator."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"outputs/top_damping_audit_20260910"
CASES={600:ROOT/"outputs/diagnostic_sequence_20260905/sequence.h5",
       300:ROOT/"outputs/resolution_300m_20260909/sequence.h5"}
CALLS_PER_STEP=4
MAX_FRACTION=.05


def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    rows=[];schedules={};case_times={}
    for dx,path in CASES.items():
        with h5py.File(path,"r") as f:
            zf=f["grid/zf"][:];nz=len(zf)-1;nd=max(2,nz//10)
            profile=np.linspace(0.,1.,nd)**2;multipliers=1-MAX_FRACTION*profile
            keys=list(f["snapshots"])
            case_times[dx]=[]
            for key in keys:
                g=f[f"snapshots/{key}"];t=float(g.attrs["time_s"])
                if not 2790<=t<=3301:continue
                w=np.asarray(g["w"][:,:,-nd:])
                before=np.sum(w*w,axis=(0,1));after=np.sum((w*multipliers[None,None,::-1])**2,axis=(0,1))
                total_before=float(before.sum());total_after=float(after.sum())
                rows.append({"grid_dx_m":dx,"time_s":t,
                             "w_damping_rms_m_s":float(np.sqrt(np.mean(w*w))),
                             "w_damping_absmax_m_s":float(np.max(np.abs(w))),
                             "instantaneous_ke_proxy_removed_fraction_per_call":
                                 float((total_before-total_after)/total_before) if total_before else 0.})
                case_times[dx].append(t)
            steps=f["steps"][:];schedules[dx]=steps
    write_csv(OUT/"snapshot_operator_effect.csv",rows)
    first=max(min(v) for v in case_times.values());last=min(max(v) for v in case_times.values())
    summary={"status":"complete","common_interval_s":[first,last],"calls_per_native_step":CALLS_PER_STEP,
             "maximum_fraction_per_call":MAX_FRACTION,
             "profile_semantics":"top face multiplier 1.0; lowest affected face multiplier 0.95",
             "cases":{},"limitations":[
                 "Applying the operator to an archived post-step state is an instantaneous counterfactual, not measured causal energy loss.",
                 "The kinetic-energy proxy omits density and staggered control-volume weights; only its within-case fraction is used.",
                 "Tendencies between boundary calls can replenish vertical velocity."]}
    for dx,steps in schedules.items():
        overlap=np.maximum(0,np.minimum(steps[:,0]+steps[:,1],last)-np.maximum(steps[:,0],first))
        active=overlap>0;calls=float(np.sum(CALLS_PER_STEP*overlap[active]/steps[active,1]))
        subset=[r for r in rows if r["grid_dx_m"]==dx]
        summary["cases"][str(dx)]={
            "overlapping_native_steps":int(active.sum()),"fractional_operator_calls":calls,
            "nominal_lowest_face_amplitude_decay_rate_s-1":float(-np.log1p(-MAX_FRACTION)*calls/(last-first)),
            "median_w_damping_rms_m_s":float(np.median([r["w_damping_rms_m_s"] for r in subset])),
            "max_w_damping_absmax_m_s":float(max(r["w_damping_absmax_m_s"] for r in subset)),
            "median_instantaneous_ke_proxy_removed_fraction_per_call":float(np.median([
                r["instantaneous_ke_proxy_removed_fraction_per_call"] for r in subset])),
        }
    a=summary["cases"]["600"];b=summary["cases"]["300"]
    summary["ratios_300_over_600"]={k:float(b[k]/a[k]) for k in (
        "nominal_lowest_face_amplitude_decay_rate_s-1","median_w_damping_rms_m_s",
        "max_w_damping_absmax_m_s","median_instantaneous_ke_proxy_removed_fraction_per_call")}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":main()
