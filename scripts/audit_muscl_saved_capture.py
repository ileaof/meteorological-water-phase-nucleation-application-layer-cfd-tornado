"""Independent scalar/HDF5 consistency audit; does not run the model."""
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = ROOT / 'outputs/muscl_direct_capture_20260920'
    out = ROOT / 'outputs/muscl_saved_audit_20261004'
    manifest = json.loads((source / 'manifest.json').read_text())
    assert sha(source / 'capture.h5') == manifest['capture_sha256']
    assert sha(source / 'metadata.json') == manifest['metadata_sha256']
    with h5py.File(source / 'capture.h5', 'r') as f:
        assert f.attrs['status'] == 'complete'
        terms = json.loads(f.attrs['terms'])
        masks = json.loads(f.attrs['masks'])
        steps = f['step_statistics'][:]
        blocks = f['block_statistics'][:]
        balances = f['block_balances'][:]
        levels = f['levels_m'][:]
        schedule = f['schedule'][:]
        for group in f['blocks'].values():
            weights = group['masks'][:]
            for mode in ('moving_end', 'moving_symmetric', 'fixed_initial'):
                ids = sorted((i for i, m in enumerate(masks) if m['mask_mode'] == mode),
                             key=lambda i: masks[i]['radius_m'])
                assert np.all(np.diff(weights[ids], axis=0) >= 0)
    errors = {}

    def check(name, a, b, atol=1e-8):
        np.testing.assert_allclose(a, b, rtol=1e-10, atol=atol)
        errors[name] = float(np.max(np.abs(a-b)))

    check('schedule_continuity_s', schedule[:-1, 0]+schedule[:-1, 1], schedule[1:, 0], 2e-8)
    assert np.all(np.diff(schedule[:, 2]) == 1)
    for prefix in ('total', 'les'):
        reconstructed = sum(steps[:, terms.index(prefix+'_flux_'+d), :, :, 0] for d in 'xyz')
        reconstructed += steps[:, terms.index(prefix+'_roundoff'), :, :, 0]
        check(prefix+'_directional_sum_m2_s', reconstructed, steps[:, terms.index(prefix+'_muscl'), :, :, 0])
    check('adjoint_step_m2_s', steps[..., 0], steps[..., 2])
    check('signed_step_block_sum_m2_s', steps[..., 0].sum(axis=0), blocks[..., 0].sum(axis=0))
    check('inventory_continuity_m2_s', balances[:-1, :, :, :, 1], balances[1:, :, :, :, 0])
    for q, names in enumerate((('total_les', 'total_muscl', 'total_other'),
                               ('les_injection', 'les_muscl'))):
        increments = sum(blocks[:, terms.index(n), :, :, 0] for n in names)
        check('inventory_balance_'+str(q)+'_m2_s', balances[:, :, q, :, 1]-balances[:, :, q, :, 0],
              increments+balances[:, :, q, :, 2])
    assert np.all(steps[..., 1]+1e-8 >= np.abs(steps[..., 0]))
    assert np.all(blocks[..., 1].sum(axis=0) <= steps[..., 1].sum(axis=0)+1e-8)
    total = steps[:, terms.index('total_muscl'), :, :, 0]
    signed = total.sum(axis=0)
    negative = (total < 0).sum(axis=0)
    rows = []
    for k, z in enumerate(levels):
        values = []
        for mode in ('moving_end', 'moving_symmetric', 'fixed_initial'):
            m = next(i for i, d in enumerate(masks) if d['radius_m'] == 4200 and d['mask_mode'] == mode)
            values.append(dict(mode=mode, signed_m2_s=float(signed[m, k]), negative_steps=int(negative[m, k])))
        rows.append(dict(level_m=float(z), modes=values))
    result = dict(status='PASS', date='2026-10-04', no_simulation=True,
                  capture_sha256=manifest['capture_sha256'], script_sha256=sha(Path(__file__)),
                  checks_max_abs=errors, nested_masks_verified=True,
                  triangle_inequalities_verified=True, levels_at_4200m=rows)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'audit.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
