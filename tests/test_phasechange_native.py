"""Native ABI, independent derivatives, linear solves and connected evolution."""
import json
from dataclasses import replace
import numpy as np
import pytest
from phasechange.finite_volume import Slab,advance
from phasechange.material import Material
from phasechange.native import NativeKernel
from phasechange.native_checks import ROOT,experiment,reference_assembly

MAT=Material(2.,3.,1.5,12.,10.,14.)


@pytest.fixture(scope='module')
def kernel():
    # A missing build is an error, never a silently skipped native test.
    return NativeKernel()


@pytest.mark.parametrize('n',[1,2,9])
@pytest.mark.parametrize('walls',[(None,None),(7.,None),(7.,19.)])
def test_independent_cell_balance_and_buffers(kernel,n,walls):
    slab=Slab(2.,n,.3); old=np.linspace(8,16,n); trial=np.linspace(16,8,n); q=np.arange(n)-2.
    saved=old.copy(); reference=reference_assembly(MAT,slab,old,trial,2.5,q,*walls)
    actual=kernel.assemble(MAT,slab,old,trial,2.5,source=q,left=walls[0],right=walls[1])
    for a,b in zip(actual,reference): np.testing.assert_allclose(a,b,rtol=2e-14,atol=2e-13)
    np.testing.assert_array_equal(old,saved)


def matrix(diagonal,off):
    return np.diag(diagonal)+np.diag(off,1)+np.diag(off,-1)


def test_centered_directional_jacobian_away_from_kinks(kernel):
    slab=Slab(2.,5,.3); old=np.full(5,12.); t=np.array([8.,11.,12.,15.,18.]); v=np.array([.3,-.8,.7,.1,-.2]); eps=1e-5
    r,d,e=kernel.assemble(MAT,slab,old,t,2.,left=7,right=19)
    plus=kernel.assemble(MAT,slab,old,t+eps*v,2.,left=7,right=19)[0]
    minus=kernel.assemble(MAT,slab,old,t-eps*v,2.,left=7,right=19)[0]
    np.testing.assert_allclose((plus-minus)/(2*eps),matrix(d,e)@v,rtol=2e-8,atol=2e-8)


@pytest.mark.parametrize('side,sign',[('right',1.),('left',-1.)])
def test_endpoint_directional_derivatives(kernel,side,sign):
    slab=Slab(2.,2,.3); t=np.array([10.,14.]); old=np.array([8.,16.]); eps=1e-5
    r,d,e=kernel.assemble(MAT,slab,old,t,2.,side=side)
    shifted=kernel.assemble(MAT,slab,old,t+sign*eps,2.)[0]
    np.testing.assert_allclose((shifted-r)/eps,matrix(d,e)@np.full(2,sign),rtol=2e-8,atol=2e-8)


@pytest.mark.parametrize('n',[1,2,17,257])
def test_ldl_manufactured_solution(kernel,n):
    d=np.linspace(2,3,n); e=np.full(n-1,-.25); expected=np.sin(np.arange(n)+.5)
    rhs=d*expected; rhs[:-1]+=e*expected[1:]; rhs[1:]+=e*expected[:-1]
    result=kernel.solve(d,e,rhs)
    np.testing.assert_allclose(result,expected,rtol=2e-14,atol=2e-14)


@pytest.mark.parametrize('d,e,rhs',[([1,1],[2],[1,1]),([0],[],[1]),([np.nan],[],[1]),([1],[1],[1]),([],[],[])])
def test_invalid_linear_system_is_rejected(kernel,d,e,rhs):
    with pytest.raises(ValueError): kernel.solve(d,e,rhs)


@pytest.mark.parametrize('args',[dict(dt=0),dict(left=-1),dict(right=np.inf),dict(side='middle'),dict(previous=-1),dict(trial=np.nan)])
def test_invalid_assembly_is_rejected(kernel,args):
    params=dict(previous=12.,trial=12.,dt=1.)
    params.update(args)
    with pytest.raises(ValueError): kernel.assemble(MAT,Slab(1.,2),**params)


def test_negative_trial_extension_and_negative_accepted_state(kernel):
    slab=Slab(1.,1)
    r,d,e=kernel.assemble(MAT,slab,[8.],[-2.],1.)
    assert r[0]==-60 and d[0]==6
    with pytest.raises(ValueError): advance(MAT,slab,8.,1.,source=-1e5,backend=kernel)


def test_native_failure_respects_iteration_limit(kernel):
    with pytest.raises(RuntimeError,match='iteration limit'):
        advance(MAT,Slab(2.,7,.3),np.linspace(8,16,7),2.,left=95.,source=10.,max_iterations=1,backend=kernel)


def test_versioned_native_evidence(kernel):
    actual=experiment(kernel); recorded=json.loads((ROOT/'results/enthalpy-native.json').read_text())
    assert actual['source_sha256']==recorded['source_sha256']
    assert len(actual['assembly_and_state_cases'])==36
    assert len(actual['connected_trajectory'])==10
    for a,b in zip(actual['connected_trajectory'],recorded['connected_trajectory']):
        np.testing.assert_allclose(a['native']['temperature'],b['native']['temperature'],rtol=0,atol=2e-9)


def test_missing_library_has_no_fallback():
    with pytest.raises(FileNotFoundError): NativeKernel(ROOT/'build/does-not-exist')
