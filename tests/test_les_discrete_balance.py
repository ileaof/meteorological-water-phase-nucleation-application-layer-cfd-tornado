"""Independent identities needed by the read-only circulation audit."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('h5py')
spec = importlib.util.spec_from_file_location('les_balance', Path(__file__).resolve().parents[1]
                                             / 'scripts/analyze_les_discrete_balance.py')
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def test_solid_rotation_circulation_and_discrete_boundary():
    x, y = np.arange(12.) * 3, np.arange(10.) * 4
    u = np.broadcast_to(-0.02*y[None,:,None], (13,10,3))
    v = np.broadcast_to(0.02*x[:,None,None], (12,11,3))
    mask = np.zeros((12,10)); mask[3:8,2:7] = 1
    zeta = analysis.curl_z(u,v,3,4)
    assert np.allclose(zeta,0.04)
    expected = np.full(3,25*12*0.04)
    assert np.allclose(analysis.integral(mask,zeta,12),expected)
    assert np.allclose(analysis.boundary_integral(mask,u,v,3,4),expected)


def test_boundary_identity_including_external_edges():
    rng = np.random.default_rng(37)
    mask = rng.random((12,10))
    u, v = rng.normal(size=(13,10,3)), rng.normal(size=(12,11,3))
    assert np.allclose(analysis.integral(mask,analysis.curl_z(u,v,3,4),12),
                       analysis.boundary_integral(mask,u,v,3,4),atol=1e-12)


@pytest.mark.parametrize('symmetric',[False,True])
def test_moving_mask_identity_and_convention(symmetric):
    rng = np.random.default_rng(17)
    qa, qb = rng.normal(size=(12,10,3)), rng.normal(size=(12,10,3))
    ma, mb = (rng.random((12,10))>.4).astype(float), (rng.random((12,10))>.4).astype(float)
    weight, motion = analysis.mask_split(qa,qb,ma,mb,symmetric)
    delta = analysis.integral(mb,qb,12)-analysis.integral(ma,qa,12)
    assert np.allclose(delta,analysis.integral(weight,qb-qa,12)+12*motion.sum(axis=(0,1)))


def test_stationary_field_change_is_entirely_mask_motion():
    q = np.arange(120.).reshape(12,10,1)
    ma = np.zeros((12,10)); ma[2:4,3:6]=1
    mb = np.zeros((12,10)); mb[3:5,3:6]=1
    _, motion = analysis.mask_split(q,q,ma,mb)
    assert np.allclose(analysis.integral(mb,q,12)-analysis.integral(ma,q,12),
                       12*motion.sum(axis=(0,1)))
