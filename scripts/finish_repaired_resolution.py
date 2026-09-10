"""Finish diagnostics after the running repaired-sequence provenance passes gates."""
from pathlib import Path
import json
import subprocess
import sys
import time
import numpy as np
import h5py

ROOT=Path(__file__).resolve().parents[1]
PROV=ROOT/'outputs/resolution_300m_provenance_20260909'
STATUS=ROOT/'outputs/resolution_300m_20260909/repaired_followup_status.json'


def record(stage,**extra):
    STATUS.write_text(json.dumps(dict(stage=stage,**extra),indent=2),encoding='utf-8')


def main():
    try:
        record('waiting_for_provenance')
        deadline=time.monotonic()+24*3600
        while not (PROV/'summary.json').exists():
            if time.monotonic()>deadline:raise TimeoutError('Provenance did not finish in 24h')
            time.sleep(30)
        # Writer may still be finishing the JSON file.
        for _ in range(10):
            try:
                s=json.loads((PROV/'summary.json').read_text(encoding='utf-8'));break
            except json.JSONDecodeError:time.sleep(1)
        assert s['status']=='complete' and s['numerical_validity_gate']=='PASS'
        assert s['native_steps']==1283 and all(r['bitwise'] for r in s['neutrality'].values())
        assert abs(s['end_time_s']-3300.023494218)<2e-8
        with h5py.File(PROV/'provenance_long_v4.h5','r') as f:
            assert f.attrs['status']=='complete'
            hist=f['closure_history'][:]
            for key in hist.dtype.names:
                if key!='stage': assert np.isfinite(hist[key]).all(),key
            limits={'omega_relative_rms':1e-10,'omega_max_abs':1e-10,
                    'tilting_relative_rms':1e-10,'tilting_max_abs':1e-12}
            for key,limit in limits.items(): assert float(hist[key].max())<=limit,key
            remainder={key:float(hist[key].max()) for key in hist.dtype.names if 'remainder' in key}
            # Conservative roundoff screening in addition to the original closure gates.
            assert max(remainder.values())<1e-10,remainder
            conditioning={key:dict(max=float(hist[key].max()),final=float(hist[key][-1]))
                          for key in ('omega_condition_index','tilting_condition_index')}
        record('numerical_gates_passed_conditioning_review_pending',conditioning=conditioning,remainder=remainder)
        subprocess.run([sys.executable,'scripts/analyze_tornadogenesis_mechanism.py',
                        '--sequence','outputs/resolution_300m_20260909/sequence.h5',
                        '--provenance',str(PROV/'provenance_long_v4.h5'),
                        '--validity',str(PROV/'summary.json'),'--radial-bin-m','600',
                        '--out','outputs/resolution_300m_audit_20260909'],cwd=ROOT,check=True)
        subprocess.run([sys.executable,'scripts/compare_resolution_audits.py'],cwd=ROOT,check=True)
        record('comparison_artifacts_complete_scientific_review_pending',conditioning=conditioning,remainder=remainder)
    except Exception as e:
        record('failed',error=repr(e));raise


if __name__=='__main__':main()
