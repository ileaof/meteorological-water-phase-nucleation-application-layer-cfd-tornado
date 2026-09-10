"""Validate the completed 300 m archive and repair only its final step attribute."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT/'outputs/resolution_300m_20260909/sequence.h5'
REPORT = PATH.parent/'final_step_repair.json'


def file_hash():
    h=hashlib.sha256()
    with PATH.open('rb') as f:
        for b in iter(lambda:f.read(16*1024**2),b''): h.update(b)
    return h.hexdigest()


def inventory(f):
    result={}
    def visit(name,obj):
        if isinstance(obj,h5py.Dataset):
            a=obj[()]
            if np.issubdtype(a.dtype,np.number) and not np.isfinite(a).all():
                raise ValueError('Nonfinite dataset '+name)
            result[name]=dict(shape=list(a.shape),dtype=str(a.dtype),
                             sha256=hashlib.sha256(a.tobytes()).hexdigest())
    f.visititems(visit)
    return result


def main():
    if REPORT.exists(): raise RuntimeError('Repair report exists; refuse repeat mutation')
    before=file_hash()
    with h5py.File(PATH,'r') as f:
        assert f.attrs['status']=='complete'
        groups=list(f['snapshots'].values()); assert len(groups)==18
        steps=f['steps'][:]; assert steps.shape==(1283,3)
        assert np.array_equal(steps[:,2],np.arange(6114,7397))
        assert (steps[:,1]>0).all()
        np.testing.assert_allclose(steps[:-1,0]+steps[:-1,1],steps[1:,0],rtol=0,atol=2e-8)
        start=float(groups[0].attrs['time_s']);end=float(groups[-1].attrs['time_s'])
        assert abs(start-steps[0,0])<2e-8
        assert abs(start+steps[:,1].sum()-end)<2e-8
        assert abs(steps[-1,0]+steps[-1,1]-end)<2e-8
        assert groups[-1].attrs['step']==7398
        times=np.r_[steps[:,0],end]
        checks=[]
        for i,g in enumerate(groups):
            assert g.attrs['complete']
            t=float(g.attrs['time_s']); j=int(np.argmin(abs(times-t)))
            assert abs(times[j]-t)<2e-8
            expected=6114+j
            if i<17: assert g.attrs['step']==expected
            else: assert expected==7397
            for n in ('u','v','w','p','p_dyn','theta','qv','ql','qi','qr','qs','qg','qh'):
                expected_shape=(240+(n=='u'),240+(n=='v'),48+(n=='w'))
                assert g[n].shape==expected_shape,(g.name,n,g[n].shape)
            if i: assert 'increments' in g and 'kinematic_integrals' in g
            checks.append(dict(snapshot=i,time_s=t,recorded_step=int(g.attrs['step']),expected_step=expected))
        datasets=inventory(f)
    report=dict(status='validated_before_repair',utc=datetime.now(timezone.utc).isoformat(),
                file_sha256_before=before,steps=1283,dt_sum_s=float(steps[:,1].sum()),
                time_start_s=start,time_end_s=end,snapshots=checks,datasets_before=datasets)
    REPORT.write_text(json.dumps(report,indent=2),encoding='utf-8')
    correction=dict(field='snapshots/00017.attrs.step',old=7398,new=7397,
                    reason='close called after driver increment; legacy save added one again',
                    file_sha256_before=before,utc=report['utc'],physical_arrays_changed=False)
    with h5py.File(PATH,'r+') as f:
        f['snapshots/00017'].attrs.modify('step',7397)
        f.attrs['metadata_corrections']=json.dumps([correction])
        f.flush()
    with h5py.File(PATH,'r') as f:
        assert inventory(f)==datasets,'Dataset integrity changed'
        assert f['snapshots/00017'].attrs['step']==7397
    report.update(status='PASS',file_sha256_after=file_hash(),all_datasets_bitwise_unchanged=True,
                  correction=correction)
    REPORT.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('datasets_before','snapshots')},indent=2))


if __name__=='__main__': main()
