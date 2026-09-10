"""Wait for the refined run, then gate provenance and generate its diagnostics."""
from pathlib import Path
import json
import subprocess
import sys
import time

import h5py

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / 'outputs/resolution_300m_20260909'
PROV = ROOT / 'outputs/resolution_300m_provenance_20260909'
AUDIT = ROOT / 'outputs/resolution_300m_audit_20260909'


def main():
    status = RUN / 'followup_status.json'
    def record(stage, **extra):
        status.write_text(json.dumps(dict(stage=stage, **extra), indent=2), encoding='utf-8')
    try:
        record('waiting_for_simulation')
        deadline = time.monotonic() + 24 * 3600
        while True:
            try:
                metadata = json.loads((RUN / 'metadata.json').read_text(encoding='utf-8'))
            except json.JSONDecodeError:
                time.sleep(5)
                continue
            if metadata['status'] == 'complete':
                break
            if metadata['status'] != 'running':
                raise RuntimeError('Refined simulation did not complete')
            if time.monotonic() > deadline:
                raise TimeoutError('Simulation did not complete within 24 hours')
            time.sleep(30)
        sequence = RUN / 'sequence.h5'
        with h5py.File(sequence, 'r') as handle:
            if handle.attrs['status'] != 'complete':
                raise RuntimeError('Sequence is not complete')
            snapshots = sorted(handle['snapshots'])
            first, last = int(snapshots[0]), int(snapshots[-1])
        record('running_provenance', first_snapshot=first, last_snapshot=last)
        subprocess.run([sys.executable, str(ROOT / 'scripts/run_vorticity_provenance_long.py'),
                        '--sequence', str(sequence), '--snapshot', str(first),
                        '--target-snapshot', str(last), '--device', 'gpu', '--out', str(PROV)],
                       cwd=ROOT, check=True)
        record('running_diagnostics')
        subprocess.run([sys.executable, str(ROOT / 'scripts/analyze_tornadogenesis_mechanism.py'),
                        '--sequence', str(sequence), '--provenance', str(PROV / 'provenance_long_v4.h5'),
                        '--validity', str(PROV / 'summary.json'), '--out', str(AUDIT)],
                       cwd=ROOT, check=True)
        record('diagnostics_complete_comparison_pending')
    except Exception as error:
        record('failed', error=repr(error))
        raise


if __name__ == '__main__':
    main()
