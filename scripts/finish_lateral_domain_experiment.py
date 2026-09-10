"""Validate and analyze the 120 km lateral-domain sequence after it finishes."""
from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RUN=ROOT/"outputs/domain_extent_120km_600m_20260910"
STATUS=RUN/"followup_status.json"


def record(stage,**extra):
    RUN.mkdir(parents=True,exist_ok=True)
    STATUS.write_text(json.dumps({"stage":stage,**extra},indent=2),encoding="utf-8")


def main():
    try:
        record("waiting_for_sequence")
        deadline=time.monotonic()+12*3600
        while True:
            metadata=RUN/"metadata.json"; sequence=RUN/"sequence.h5"
            if metadata.exists() and sequence.exists():
                m=json.loads(metadata.read_text(encoding="utf-8"))
                if m.get("status")=="complete": break
                if m.get("status")=="interrupted": raise RuntimeError("simulation interrupted")
            if time.monotonic()>deadline: raise TimeoutError("simulation did not finish in 12 h")
            time.sleep(30)
        with h5py.File(sequence,"r") as f:
            assert f.attrs["status"]=="complete"
            assert (len(f["grid/xc"]),len(f["grid/yc"]),len(f["grid/zc"]))==(200,200,48)
            assert abs(float(f["grid/xf"][-1])-120000.)<1e-9
            assert abs(float(f["grid/yf"][-1])-120000.)<1e-9
            assert abs(float(f["grid/zf"][-1])-15000.)<1e-9
            keys=list(f["snapshots"]); assert len(keys)>=18,len(keys)
            times=np.array([float(f[f"snapshots/{k}"].attrs["time_s"]) for k in keys])
            steps=np.array([int(f[f"snapshots/{k}"].attrs["step"]) for k in keys])
            assert np.all(np.diff(times)>0) and np.all(np.diff(steps)>0)
            assert abs(times[-1]-3300.023494218)<2e-8
            for key in keys:
                g=f[f"snapshots/{key}"]
                for name in ("u","v","w","p","p_dyn","theta","qv","ql","qi","qr","qs","qg","qh","rho"):
                    assert np.isfinite(g[name][:]).all(),(key,name)
        record("sequence_validated",frames=len(keys),start_time_s=float(times[0]),end_time_s=float(times[-1]))
        subprocess.run([sys.executable,"scripts/audit_domain_extent.py","--sequence",str(sequence),
                        "--out","outputs/domain_extent_120km_audit_20260910"],cwd=ROOT,check=True)
        subprocess.run([sys.executable,"scripts/compare_lateral_domain_sequences.py"],cwd=ROOT,check=True)
        result=json.loads((ROOT/"outputs/domain_extent_120km_comparison_20260910/summary.json").read_text())
        record("comparison_complete_scientific_review_pending",frames=len(keys),
               screening_classification=result["screening_classification"],
               materially_favorable_primary_metrics=result["materially_favorable_primary_metrics"])
    except Exception as exc:
        record("failed",error=repr(exc));raise


if __name__=="__main__": main()
