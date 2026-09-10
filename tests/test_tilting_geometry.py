"""Analytic checks of the offline geometric attribution, without a simulation."""
import importlib.util
from pathlib import Path
import numpy as np
import pytest

pytest.importorskip('h5py')
spec=importlib.util.spec_from_file_location('tilting',Path(__file__).resolve().parents[1]/'scripts/analyze_tilting_reversal.py')
analysis=importlib.util.module_from_spec(spec); spec.loader.exec_module(analysis)


def test_exact_geometric_attribution_arbitrary_vectors():
    rng=np.random.default_rng(18)
    for _ in range(100):
        vectors=rng.normal(size=(4,2))
        result=analysis.angular_decomposition(*vectors)
        assert abs(result['closure'])<2e-14


@pytest.mark.parametrize('rotating',['omega','gradient'])
def test_only_one_vector_rotates(rotating):
    start=np.array([1.,0.]); end=np.array([-.5,np.sqrt(.75)])
    vectors=(start,end,start,start) if rotating=='omega' else (start,start,start,end)
    result=analysis.angular_decomposition(*vectors)
    assert result[rotating+'_orientation']==pytest.approx(-1.5)
    other='gradient' if rotating=='omega' else 'omega'
    assert result[other+'_orientation']==pytest.approx(0.)
    assert result['omega_magnitude']==pytest.approx(0.)
    assert result['gradient_magnitude']==pytest.approx(0.)


def test_trilinear_derivative_nonuniform_coordinates():
    coords=(np.array([0.,2.,5.]),np.array([-2.,1.,7.]),np.array([0.,.3,1.7,4.]))
    x,y,z=np.meshgrid(*coords,indexing='ij')
    values=2*x-3*y+4*z+.7*x*y*z
    points=np.array([[1.,0.,.2],[3.,4.,2.]])
    xp,yp,zp=points.T
    expected=np.column_stack([2+.7*yp*zp,-3+.7*xp*zp,4+.7*xp*yp])
    np.testing.assert_allclose(analysis.trilinear_gradient(values,coords,points),expected,atol=1e-12)
