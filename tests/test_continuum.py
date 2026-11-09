"""Independent PDE, interface-flux and matrix-exponential reference checks."""
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pytest
from scipy.linalg import expm

from phasechange.continuum import similarity_parameter, stefan_solution, cosine_reference
from phasechange.material import Material
from phasechange.refinement import heat_study, acceptance

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('ste',[1e-300,1e-20,1e20,1e300])
def test_similarity_extreme_scales_in_log_equation(ste):
    lam = similarity_parameter(ste)
    assert math.log(lam)+math.log(math.erf(lam))+lam*lam == pytest.approx(
        math.log(ste)-.5*math.log(math.pi),abs=2e-12)


@pytest.mark.parametrize('lam',[.01,.1,.5,1.,2.])
def test_similarity_root_for_prescribed_interface_constant(lam):
    ste = math.sqrt(math.pi)*lam*math.exp(lam*lam)*math.erf(lam)
    assert similarity_parameter(ste) == pytest.approx(lam,abs=1e-14)


@pytest.mark.parametrize('ste',[0.,-1.,np.inf,np.nan])
def test_invalid_similarity_number(ste):
    with pytest.raises(ValueError):
        similarity_parameter(ste)


def test_reference_heat_equation_wall_and_stefan_jump():
    # Independent differential identities, using point differences of T.
    rho,c,k,latent,tm,wall = 1200.,900.,1.8,90000.,310.,330.
    alpha = k/(rho*c)
    kw = dict(diffusivity=alpha,stefan_number=c*(wall-tm)/latent,melting=tm,wall=wall)
    t = 60.
    front,_ = stefan_solution([],t,**kw)
    x = .4*front; dx = front*1e-3; dt = t*1e-3
    _,temps = stefan_solution([x-dx,x,x+dx],t,**kw)
    _,early = stefan_solution([x],t-dt,**kw)
    _,late = stefan_solution([x],t+dt,**kw)
    assert (late[0]-early[0])/(2*dt) == pytest.approx(
        alpha*(temps[0]-2*temps[1]+temps[2])/dx**2,rel=3e-6)
    _,bounds = stefan_solution([0,front,2*front],t,**kw)
    np.testing.assert_allclose(bounds,[wall,tm,tm],atol=1e-12,rtol=0)
    _,inside = stefan_solution([front-dx,front-2*dx],t,**kw)
    left_derivative = (3*tm-4*inside[0]+inside[1])/(2*dx)
    assert -k*left_derivative == pytest.approx(rho*latent*front/(2*t),rel=2e-6)
    later,_ = stefan_solution([],4*t,**kw)
    assert later == pytest.approx(2*front,rel=1e-14)


@pytest.mark.parametrize('kwargs',[dict(time=0),dict(time=-1),dict(diffusivity=0),
    dict(wall=310),dict(x=[-1]),dict(x=[np.nan])])
def test_invalid_stefan_state(kwargs):
    args = dict(x=[0],time=1,diffusivity=1,stefan_number=.2,melting=310,wall=330)
    args.update(kwargs)
    with pytest.raises(ValueError):
        stefan_solution(**args)


@pytest.mark.parametrize('cells',[2,5,17])
def test_semidiscrete_cosine_against_independent_matrix_exponential(cells):
    length,alpha,time = 2.,.3,.7
    # Sum incidence outer products, without the oracle's eigenvalue formula.
    incidence = np.zeros((cells-1,cells))
    for i in range(cells-1):
        incidence[i,i:i+2] = [-1,1]
    operator = -alpha/(length/cells)**2*(incidence.T@incidence)
    initial = cosine_reference(cells,length,alpha,0)
    actual = expm(time*operator)@initial
    np.testing.assert_allclose(actual,cosine_reference(cells,length,alpha,time,discrete=True),atol=2e-14,rtol=0)
    edges = np.linspace(0,length,cells+1)
    integrated = np.diff(np.sin(np.pi*edges/length))/(np.pi/cells)
    np.testing.assert_allclose(initial,integrated,atol=2e-15,rtol=0)


@pytest.mark.parametrize('native',[False,True])
def test_zero_latent_heat_isolated_orders(native):
    from phasechange.native import NativeKernel
    material = Material(**json.loads((ROOT/'phasechange/benchmarks.json').read_text())['material'])
    data = heat_study(material,NativeKernel() if native else None)
    assert all(.9<=r['order']<=1.1 for r in data['temporal'][1:])
    assert all(1.8<=r['order']<=2.2 for r in data['spatial'][1:])
    assert max(r['algebraic_error_K'] for r in data['temporal'])<1e-9


def test_fixed_protocol_evidence_and_failed_front_gate():
    data = json.loads((ROOT/'results/enthalpy-refinement.json').read_text())
    protocol = json.loads((ROOT/'phasechange/benchmarks.json').read_text())['stefan_plan']
    assert data['protocol']==protocol
    assert data['acceptance']==acceptance(data)
    for path,digest in data['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    for case in data['cases'].values():
        assert case['enthalpy']['cells']==case['frozen']['cells']
        for result in case.values():
            assert [s['time_s'] for s in result['samples']]==protocol['observation_times_s']
    for bad in (None,.021):
        data['cases']['refined']['enthalpy']['samples'][0]['relative_front_error']=bad
        assert not acceptance(data)['refined_stefan_front']
