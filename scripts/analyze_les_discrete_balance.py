"""Read-only, coincident-interval LES-label and total circulation accounting.

Definitions and units: docs/LES_DISCRETE_BALANCE_METHOD.md.
This module never imports or advances a simulation.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUNS = {
    600: ('diagnostic_sequence_20260905', 'vorticity_provenance_long_v4_2790_3300',
          'resolution_600m_comparison_audit_20260909'),
    300: ('resolution_300m_20260909', 'resolution_300m_provenance_20260909',
          'resolution_300m_audit_20260909'),
}
ENDPOINTS = (0, 1, 3, 4, 8, 9, 11, 14, 15, 17)
RADII = (1200.0, 2400.0, 4200.0, 6000.0)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tracer_version(path, recorded_sha):
    if sha(path) == recorded_sha:
        return 'exact current file'
    # Reconstruct only the documented close-time step-counter correction in memory.
    lines = path.read_bytes().splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines) if b'is called after the outer driver' in line)
    end = next(i for i in range(start, len(lines)) if b'self.save(sim, initial=True)' in lines[i])
    old = (b''.join(lines[:start])
           + lines[end].replace(b'self.save(sim, initial=True)', b'self.save(sim)')
           + b''.join(lines[end+1:]))
    if hashlib.sha256(old).hexdigest() != recorded_sha:
        raise ValueError('unverified tracer version; cannot identify LES-label evolution')
    return 'exact recorded SHA after in-memory reversal of close-time step-counter fix only'


def write_csv(path, rows):
    with path.open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class Archive:
    def __init__(self, path):
        self.path = path
        self.file = h5py.File(path, 'r')
        self.digest = hashlib.sha256()
        self.selections = []
        self.before = path.stat()

    def read(self, name, selection=...):
        array = np.asarray(self.file[name][selection])
        descriptor = dict(dataset=name, selection=repr(selection),
                          shape=list(array.shape), dtype=array.dtype.str)
        self.digest.update(json.dumps(descriptor, sort_keys=True).encode())
        self.digest.update(np.ascontiguousarray(array).tobytes())
        self.selections.append(descriptor)
        if not np.isfinite(array).all():
            raise ValueError(f'nonfinite input: {self.path}:{name}')
        return array

    def close(self):
        metadata = {k: str(v) for k, v in self.file.attrs.items()}
        self.file.close()
        after = self.path.stat()
        if (after.st_size, after.st_mtime_ns) != (self.before.st_size, self.before.st_mtime_ns):
            raise RuntimeError(f'input changed during analysis: {self.path}')
        return dict(path=str(self.path.relative_to(ROOT)), bytes=after.st_size,
                    mtime_ns=after.st_mtime_ns, attributes=metadata,
                    selected_data_sha256=self.digest.hexdigest(),
                    hash_semantics='ordered descriptor JSON followed by C-order selected array bytes',
                    selections=self.selections)


def derivative(a, spacing, axis):
    return np.gradient(a, spacing, axis=axis, edge_order=1)


def derivative_transpose(mask, spacing, axis):
    m = np.moveaxis(mask, axis, 0)
    result = np.zeros_like(m, dtype=float)
    result[0] -= m[0] / spacing
    result[1] += m[0] / spacing
    result[-2] -= m[-1] / spacing
    result[-1] += m[-1] / spacing
    result[:-2] -= m[1:-1] / (2 * spacing)
    result[2:] += m[1:-1] / (2 * spacing)
    return np.moveaxis(result, 0, axis)


def center_uv(u, v):
    return (u[:-1] + u[1:]) * 0.5, (v[:, :-1] + v[:, 1:]) * 0.5


def curl_z(u, v, dx, dy):
    uc, vc = center_uv(u, v)
    return derivative(vc, dx, 0) - derivative(uc, dy, 1)


def integral(mask, field, area):
    return area * np.einsum('ij,ijk->k', mask, field, optimize=True)


def boundary_integral(mask, u, v, dx, dy):
    uc, vc = center_uv(u, v)
    return (integral(derivative_transpose(mask, dx, 0), vc, dx * dy)
            - integral(derivative_transpose(mask, dy, 1), uc, dx * dy))


def mask_split(qa, qb, ma, mb, symmetric=False):
    weight = 0.5 * (ma + mb) if symmetric else mb
    qmask = 0.5 * (qa + qb) if symmetric else qa
    return weight, (mb - ma)[:, :, None] * qmask


def exact_matches(seq, prov):
    skeys = sorted(seq['snapshots'])
    stimes = np.array([seq[f'snapshots/{k}'].attrs['time_s'] for k in skeys])
    result = {}
    for pkey in sorted(prov['snapshots']):
        pg = prov[f'snapshots/{pkey}']
        matches = np.flatnonzero(np.abs(stimes - float(pg.attrs['time_s'])) <= 1e-8)
        if len(matches) == 1:
            si = int(matches[0])
            sg = seq[f'snapshots/{skeys[si]}']
            if int(sg.attrs['step']) != int(pg.attrs['step']):
                raise ValueError('time coincidence without step coincidence')
            result[int(pkey)] = si
    return skeys, result


def summarize_terms(rows):
    groups = defaultdict(list)
    keys = ('grid_dx_m', 'radius_m', 'mask_mode', 'level_m', 'quantity', 'term')
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    result = []
    for key, group in sorted(groups.items()):
        signed = sum(r['signed_m2_s'] for r in group)
        temporal_abs = sum(abs(r['signed_m2_s']) for r in group)
        space_abs = sum(r['spatial_abs_m2_s'] for r in group)
        result.append(dict(zip(keys, key), blocks=len(group), signed_m2_s=signed,
                           sum_abs_block_integrals_m2_s=temporal_abs,
                           sum_spatial_abs_blocks_m2_s=space_abs,
                           spatial_retained_fraction=temporal_abs / space_abs if space_abs else None,
                           temporal_retained_fraction=abs(signed) / temporal_abs if temporal_abs else None))
    return result


def analyze_case(dx, terms, balances, inputs):
    sequence, provenance, audit = RUNS[dx]
    seq = Archive(ROOT / 'outputs' / sequence / 'sequence.h5')
    prov = Archive(ROOT / 'outputs' / provenance / 'provenance_long_v4.h5')
    track_path = ROOT / 'outputs' / audit / 'vortex_timeseries.csv'
    with track_path.open(encoding='utf-8') as stream:
        track = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(stream)]
    metadata_path = ROOT / 'outputs' / provenance / 'summary.json'
    metadata = json.loads(metadata_path.read_text())
    assert metadata['numerical_validity_gate'] == 'PASS'
    assert all(v['bitwise'] for v in metadata['neutrality'].values())
    tracer_path = ROOT / 'src/storm_dynamics/vorticity_provenance.py'
    stored_tracer_sha = next(v for k, v in metadata['metadata']['source_sha256'].items()
                             if k.replace('\\', '/').endswith('/vorticity_provenance.py'))
    version_verification = tracer_version(tracer_path, stored_tracer_sha)
    nk = int(prov.file.attrs['analysis_nz'])
    x, y, z = [seq.read('grid/' + n) for n in ('xc', 'yc', 'zc')]
    for n, expected in (('xc', x), ('yc', y), ('zc', z[:int(prov.file.attrs['output_nz_with_halo'])])):
        assert np.array_equal(prov.read('grid/' + n), expected)
    sx, sy = float(x[1] - x[0]), float(y[1] - y[0])
    assert sx == sy == dx
    skeys, matches = exact_matches(seq.file, prov.file)
    assert all(e in matches for e in ENDPOINTS)
    schedule = seq.read('steps')
    max_state_error = max_curl_closure = max_boundary_error = max_identity_error = 0.0
    block_rows = []

    def state(endpoint):
        pkey = f'snapshots/{endpoint:05d}'
        skey = f'snapshots/{skeys[matches[endpoint]]}'
        t = float(prov.file[pkey].attrs['time_s'])
        row = min(track, key=lambda r: abs(r['time_s'] - t))
        assert abs(row['time_s'] - t) < 1e-8
        les = prov.read(pkey + '/omega_source/les', (2, slice(None), slice(None), slice(0, nk)))
        total = prov.read(pkey + '/omega_total', (2, slice(None), slice(None), slice(0, nk)))
        u, v = [seq.read(skey + '/' + n, (slice(None), slice(None), slice(0, nk))) for n in ('u', 'v')]
        error = float(np.max(np.abs(curl_z(u, v, sx, sy) - total)))
        if error > 1e-12:
            raise ValueError(f'archive/replay zeta mismatch {error}')
        return dict(time=t, les=les, total=total, center=(row['center_x_m'], row['center_y_m']),
                    step=int(seq.file[skey].attrs['step']), error=error)

    def disc(center, radius):
        mask = (((x[:, None] - center[0])**2 + (y[None, :] - center[1])**2) <= radius**2).astype(float)
        assert not (mask[:2].any() or mask[-2:].any() or mask[:, :2].any() or mask[:, -2:].any())
        return mask

    a = state(ENDPOINTS[0])
    origin = a['center']
    max_state_error = a['error']
    for block, (pa, pb) in enumerate(zip(ENDPOINTS, ENDPOINTS[1:])):
        b = state(pb)
        max_state_error = max(max_state_error, b['error'])
        native = schedule[(schedule[:, 2] >= a['step']) & (schedule[:, 2] < b['step'])]
        assert len(native) == b['step'] - a['step']
        assert np.array_equal(native[:, 2], np.arange(a['step'], b['step']))
        assert abs(native[0, 0] - a['time']) < 1e-8
        assert abs(native[-1, 0] + native[-1, 1] - b['time']) < 1e-8
        assert np.max(abs(native[:-1, 0] + native[:-1, 1] - native[1:, 0])) < 1e-8
        inc = {}
        expected_start = a['time']
        for si in range(matches[pa] + 1, matches[pb] + 1):
            base = 'snapshots/' + skeys[si]
            sg = seq.file[base]
            assert abs(float(sg.attrs['interval_start_s']) - expected_start) < 1e-8
            expected_start = float(sg.attrs['time_s'])
            for op in sg['increments']:
                uv = [seq.read(base + '/increments/' + op + '/' + n,
                               (slice(None), slice(None), slice(0, nk))) for n in ('u', 'v')]
                if op not in inc:
                    inc[op] = uv
                else:
                    for accum, value in zip(inc[op], uv):
                        accum += value
        assert abs(expected_start - b['time']) < 1e-8
        fields = {op: curl_z(*uv, sx, sy) for op, uv in inc.items()}
        total_delta = b['total'] - a['total']
        closure = total_delta - sum(fields.values())
        rel = float(np.linalg.norm(closure.ravel()) / max(np.linalg.norm(total_delta.ravel()), 1e-30))
        max_curl_closure = max(max_curl_closure, rel)
        assert rel < 1e-10
        label_evolution = b['les'] - a['les'] - fields['les']
        block_rows.append(dict(block=block, provenance_start=pa, provenance_end=pb,
                               time_start_s=a['time'], time_end_s=b['time'],
                               duration_s=b['time']-a['time'], native_steps=len(native),
                               source_intervals=matches[pb]-matches[pa],
                               total_curl_closure_relative_rms=rel))
        for radius in RADII:
            moving_a, moving_b = disc(a['center'], radius), disc(b['center'], radius)
            fixed = disc(origin, radius)
            for mode in ('moving_end', 'moving_symmetric', 'fixed_initial'):
                ma, mb = (fixed, fixed) if mode == 'fixed_initial' else (moving_a, moving_b)
                symmetric = mode == 'moving_symmetric'
                weight = (ma + mb) * 0.5 if symmetric else mb
                op_integrals = {op: integral(weight, value, sx*sy) for op, value in fields.items()}
                boundary = {op: boundary_integral(weight, *uv, sx, sy) for op, uv in inc.items()}
                for op in fields:
                    err = np.max(np.abs(op_integrals[op] - boundary[op]))
                    max_boundary_error = max(max_boundary_error, float(err))
                    assert np.allclose(op_integrals[op], boundary[op], atol=1e-7, rtol=1e-10)
                for quantity in ('les_label', 'total'):
                    qkey = 'les' if quantity == 'les_label' else 'total'
                    qa, qb = a[qkey], b[qkey]
                    _, motion_field = mask_split(qa, qb, ma, mb, symmetric)
                    term_fields = ({'local_les': fields['les'], 'inferred_muscl_label': label_evolution}
                                   if quantity == 'les_label' else fields)
                    measured = {}
                    for term, field in term_fields.items():
                        signed = integral(weight, field, sx*sy)
                        absolute = integral(weight, np.abs(field), sx*sy)
                        boundary_value = boundary.get('les' if term == 'local_les' else term)
                        measured[term] = signed
                        for k in range(nk):
                            terms.append(dict(grid_dx_m=dx, block=block, radius_m=radius,
                                mask_mode=mode, level_m=float(z[k]), quantity=quantity, term=term,
                                signed_m2_s=float(signed[k]), spatial_abs_m2_s=float(absolute[k]),
                                boundary_curl_integral_m2_s=(float(boundary_value[k]) if boundary_value is not None else None)))
                    motion = sx * sy * motion_field.sum(axis=(0, 1))
                    motion_abs = sx * sy * np.abs(motion_field).sum(axis=(0, 1))
                    ia, ib = integral(ma, qa, sx*sy), integral(mb, qb, sx*sy)
                    error = ib - ia - sum(measured.values()) - motion
                    max_identity_error = max(max_identity_error, float(np.max(abs(error))))
                    assert np.allclose(ib-ia, sum(measured.values())+motion, atol=1e-7, rtol=1e-10)
                    for k in range(nk):
                        common = dict(grid_dx_m=dx, block=block, radius_m=radius, mask_mode=mode,
                                      level_m=float(z[k]), quantity=quantity)
                        terms.append(dict(**common, term='mask_motion', signed_m2_s=float(motion[k]),
                            spatial_abs_m2_s=float(motion_abs[k]), boundary_curl_integral_m2_s=None))
                        balances.append(dict(**common, time_start_s=a['time'], time_end_s=b['time'],
                            inventory_start_m2_s=float(ia[k]), inventory_end_m2_s=float(ib[k]),
                            delta_inventory_m2_s=float(ib[k]-ia[k]),
                            fixed_mask_evolution_m2_s=float(sum(measured.values())[k]),
                            mask_motion_m2_s=float(motion[k]), closure_error_m2_s=float(error[k])))
        print(f'{dx} m: block {block+1}/9, {a["time"]:.3f}-{b["time"]:.3f} s, closure={rel:.2e}', flush=True)
        a = b
    inputs.extend([seq.close(), prov.close()])
    context = {}
    for metric in ('zeta_max_s-1', 'circulation_1200_m2_s', 'circulation_4200_m2_s',
                   'vtheta_max_m_s', 'zeta_halfmax_width_m', 'rmw_cells',
                   'axis_tilt_surface_to_2km_m', 'axis_tilt_max_m',
                   'convergence_max_s-1', 'corr_zeta_convergence'):
        values = np.array([r[metric] for r in track])
        context[metric] = dict(initial=float(values[0]), final=float(values[-1]),
            min=float(values.min()), max=float(values.max()), median=float(np.median(values)),
            maximum_time_s=track[int(np.argmax(values))]['time_s'])
    return dict(exact_endpoint_indices=sorted(matches), used_endpoints=list(ENDPOINTS), blocks=block_rows,
                archive_replay_max_zeta_difference_s_1=max_state_error,
                total_curl_closure_max_relative_rms=max_curl_closure,
                boundary_identity_max_abs_m2_s=max_boundary_error,
                moving_balance_max_abs_error_m2_s=max_identity_error,
                context_from_existing_tracking=context, track_sha256=sha(track_path),
                archived_gate_summary_sha256=sha(metadata_path),
                tracer_source_version_verification=version_verification,
                archived_provenance_full_file_sha256_not_recomputed=metadata['archive']['sha256'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT/'outputs/les_discrete_balance_20260915_v2')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    terms, balances, inputs, cases = [], [], [], {}
    for dx in (600, 300):
        cases[str(dx)] = analyze_case(dx, terms, balances, inputs)
    write_csv(args.out/'terms_by_block.csv', terms)
    write_csv(args.out/'balances_by_block.csv', balances)
    write_csv(args.out/'terms_by_height.csv', summarize_terms(terms))
    cross_time = max(abs(a[k]-b[k]) for a,b in zip(cases['600']['blocks'],cases['300']['blocks'])
                     for k in ('time_start_s','time_end_s'))
    summary = dict(status='complete', created_utc=datetime.now(timezone.utc).isoformat(),
        analysis='postprocessing only; no simulation', cases=cases,
        max_cross_resolution_endpoint_offset_s=cross_time,
        units='horizontal circulation and accumulated increments: m2/s; local fields: 1/s',
        mask_geometry='same low-level center at every height; moving end, symmetric, or fixed initial',
        les_evolution_semantics='inferred by exact label difference minus local injection; not independent flux',
        boundary_semantics='transpose of diagnostic curl applied to mask, not advective zeta flux',
        wall_clock_s=time.perf_counter()-start)
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    source_paths = [Path(__file__), ROOT/'docs/LES_DISCRETE_BALANCE_METHOD.md',
                    ROOT/'src/storm_dynamics/vorticity_provenance.py', ROOT/'src/meteorological_flow/grid.py',
                    ROOT/'src/storm_dynamics/diagnostic_capture.py']
    manifest = dict(inputs=inputs, source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_paths},
                    git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                    artifacts={p.name:dict(bytes=p.stat().st_size,sha256=sha(p)) for p in args.out.iterdir() if p.is_file()})
    (args.out/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(f'Complete: {args.out}, {summary["wall_clock_s"]:.1f} s', flush=True)


if __name__ == '__main__':
    main()
