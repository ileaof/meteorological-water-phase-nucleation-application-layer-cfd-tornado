"""Read-only geometric screening of the archived 600 m storm; no causal claim."""
from pathlib import Path
import csv
import argparse
import json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequence',type=Path,default=ROOT/'outputs/diagnostic_sequence_20260905/sequence.h5')
    parser.add_argument('--out',type=Path,default=ROOT/'outputs/domain_extent_audit')
    args=parser.parse_args();source=args.sequence;out=args.out
    out.mkdir(exist_ok=True)
    rows = []
    with h5py.File(source, 'r') as f:
        x, y, z, zf = (f['grid/' + n][:] for n in ('xc', 'yc', 'zc', 'zf'))
        xf, yf = f['grid/xf'][:], f['grid/yf'][:]
        cell_area=float((x[1]-x[0])*(y[1]-y[0]))
        edge = np.minimum(np.minimum(x-xf[0], xf[-1]-x)[:, None],
                          np.minimum(y-yf[0], yf[-1]-y)[None, :])
        base_tv = f['base/theta0'][:] * (1 + .61*f['base/qv0'][:])
        k = int(np.argmin(abs(z-100)))
        nd = max(2, len(z)//10)
        damp_bottom = float(zf[-nd])
        for name, g in f['snapshots'].items():
            t = float(g.attrs['time_s'])
            if not 2790 <= t <= 3301:
                continue
            wface = g['w'][:]
            w = (wface[:, :, :-1]+wface[:, :, 1:])/2
            condensate = sum(g[n][:] for n in ('ql','qi','qr','qs','qg','qh'))
            tv = g['theta'][:,:,k]*(1+.61*g['qv'][:,:,k]-condensate[:,:,k])
            cold = tv-base_tv[k] < -1
            row = dict(snapshot=name, time_s=t, w_max_ms=float(w.max()),
                       w_top_face_absmax_ms=float(abs(wface[:,:,-1]).max()),
                       w_above_12km_max_ms=float(w[:,:,z>=12000].max()),
                       w_damping_region_max_ms=float(w[:,:,z>=damp_bottom].max()))
            masks = {'cold_m1K': cold, 'updraft_5ms': (w>5).any(axis=2),
                     'updraft_10ms': (w>10).any(axis=2)}
            for threshold in (1e-5,1e-4,1e-3):
                cloud = condensate>threshold
                key = f'condensate_{threshold:g}'
                masks[key] = cloud.any(axis=2)
                occupied = cloud.any(axis=(0,1))
                row[key+'_top_m'] = float(z[occupied].max()) if occupied.any() else None
                row[key+'_top_level_area_km2'] = float(cloud[:,:,-1].sum()*cell_area/1e6)
            for key, mask in masks.items():
                row[key+'_edge_distance_m'] = float(edge[mask].min()) if mask.any() else None
                row[key+'_area_km2'] = float(mask.sum()*cell_area/1e6)
                row[key+'_area_within_6km_edge_km2'] = float((mask & (edge<6000)).sum()*cell_area/1e6)
            rows.append(row)
            print(f'audited {name} t={t:.1f}', flush=True)
    with (out/'timeseries.csv').open('w', newline='', encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    ranges={}
    for key in rows[0]:
        if key=='snapshot': continue
        values=[r[key] for r in rows if r[key] is not None]
        ranges[key]={'min':min(values),'max':max(values)} if values else None
    summary=dict(source=str(source), frames=len(rows), domain_km=[float(xf[-1]-xf[0])/1000,float(yf[-1]-yf[0])/1000,float(zf[-1]-zf[0])/1000],
                 damping_lowest_affected_face_m=damp_bottom,
                 top_cell_center_m=float(z[-1]), ranges=ranges,
                 limitations=['Geometric screening only; no domain-sensitivity experiment.',
                               'All structures included, not only the tracked supercell.',
                               'Periodic seam contact does not prove self-interaction.',
                               'Condensate includes all six cloud/precipitation species.'])
    (out/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')


if __name__=='__main__':
    main()
