"""Independent frozen-step equations, misleading energy cancellation, and fronts."""
from dataclasses import replace
import hashlib
import json
import subprocess
import sys
import numpy as np
import pytest

from phasechange.comparison import ROOT, experiment, front_position, trajectory, verify_reproduction
from phasechange.finite_volume import Slab, advance
from phasechange.frozen import advance_frozen
from phasechange.material import Material
from phasechange.native import NativeKernel

MAT = Material(2., 3., 1.5, 12., 10., 14.)


@pytest.mark.parametrize('initial,power,dt,expected,defect', [
    (8., 12., 10., 28., 24.), (24., -12., 10., 4., -24.),
    (11., 6., 1., 11.5, 0.), (8., 6., 1., 9., 0.),
    (10., -12., 1., 9., 6.), (14., -12., 1., 12., -12.)])
@pytest.mark.parametrize('native', [False, True])
def test_exact_frozen_update_and_true_latent_defect(initial, power, dt, expected, defect, native):
    step = advance_frozen(MAT, Slab(1., 1), initial, dt, source=power,
                          backend=NativeKernel() if native else None)
    assert step.temperature == pytest.approx((expected,))
    assert step.ledger.defect_J == pytest.approx(defect, abs=1e-12)
    assert step.scaled_residual < 1e-12


@pytest.mark.parametrize('native', [False, True])
def test_independent_two_cell_mixed_wall_and_source_equation(native):
    # dx=A=dt=1; capacities 6 and 12, face=1.5, left half-face=3.
    # Absolute-temperature balance: [[10.5,-1.5],[-1.5,13.5]] T = [80,140].
    initial = np.array([8., 12.]); copy = initial.copy()
    step = advance_frozen(MAT, Slab(2., 2), initial, 1., source=[2., -4.], left=10.,
                          backend=NativeKernel() if native else None)
    determinant = 10.5*13.5-1.5**2
    expected = [(80*13.5+1.5*140)/determinant, (10.5*140+1.5*80)/determinant]
    np.testing.assert_allclose(step.temperature, expected, rtol=0, atol=1e-12)
    np.testing.assert_array_equal(initial, copy)


@pytest.mark.parametrize('cells', [1, 2, 19])
@pytest.mark.parametrize('native', [False, True])
def test_zero_latent_heat_methods_tie(cells, native):
    mat = replace(MAT, latent_heat=0.)
    mesh = Slab(2., cells)
    kwargs = dict(source=3., left=20., right=8., backend=NativeKernel() if native else None)
    initial = np.linspace(9., 16., cells)
    actual = advance_frozen(mat, mesh, initial, .3, **kwargs)
    reference = advance(mat, mesh, initial, .3, **kwargs)
    np.testing.assert_allclose(actual.temperature, reference.temperature, atol=1e-12, rtol=0)
    assert actual.ledger.scaled_defect < 1e-12


@pytest.mark.parametrize('kwargs', [dict(dt=0), dict(dt=np.inf), dict(left=-1),
    dict(right=np.nan), dict(source=[1,2]), dict(previous=-1), dict(source=-1e5)])
def test_invalid_and_negative_kelvin_results_rejected(kwargs):
    args = dict(previous=12., dt=1.)
    args.update(kwargs)
    with pytest.raises(ValueError):
        advance_frozen(MAT, Slab(1., 1), **args)


def test_closed_heat_cycle_does_not_hide_transient_defect():
    data = trajectory(MAT, Slab(1., 1), 8., [(10.,None,None,12.),(10.,None,None,-12.)], method='frozen')
    assert data['samples'][-1]['temperature'] == (8.,)
    assert data['samples'][-1]['cumulative_defect_J'] == pytest.approx(0.)
    assert data['max_scaled_cumulative_energy_defect'] == pytest.approx(.5)
    assert data['max_scaled_ledger_closure'] < 1e-12


@pytest.mark.parametrize('fraction,position,status', [
    ([1.,.8,.2,0.], 2., 'resolved'), ([1.,.5,.1,0.], 1.5, 'resolved'),
    ([1.,.5,.5,0.], None, 'plateau'), ([0.,0.,0.,0.], None, 'outside_cell_centres'),
    ([1.,1.,1.,1.], None, 'outside_cell_centres'), ([1.,.3,.7,0.], None, 'nonmonotone')])
def test_front_definition_and_unresolved_crossings(fraction, position, status):
    actual = front_position(fraction, Slab(4.,4))
    assert actual == {'position_m': position, 'status': status}


@pytest.mark.parametrize('fraction', [[1.,-1.], [1.,np.nan], [1.]])
def test_invalid_front_fields(fraction):
    with pytest.raises(ValueError):
        front_position(fraction, Slab(2.,2))


def test_complete_source_bound_experiment_reproduces():
    recorded = json.loads((ROOT/'results/enthalpy-comparison.json').read_text())
    for path, digest in recorded['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == digest
    actual = experiment()
    verify_reproduction(actual, recorded)
    assert (len(actual['uniform']),len(actual['heat_cool_slab']),len(actual['front_comparison'])) == (30,5,4)
    ties = [r for r in actual['uniform'] if r['case'] == 'solid_heating']
    assert all(r['methods']['frozen']['max_temperature_error_K'] < 1e-10 for r in ties)
    crossing = next(r for r in actual['uniform'] if r['case']=='cross_both_heating' and r['subdivisions']==1)
    assert crossing['methods']['frozen']['max_temperature_error_K'] == pytest.approx(100.)


def test_cli_rejects_native_baseline_instead_of_silent_fallback():
    result = subprocess.run([sys.executable,'-m','phasechange','--backend','native'],
                            cwd=ROOT, capture_output=True, text=True)
    assert result.returncode == 2
    assert 'choose --study comparison' in result.stderr


@pytest.mark.parametrize('key,value', [('temperature',[310.01]), ('liquid_fraction',[.51]),
    ('position_m',.002), ('stored_J',100.), ('temperature',[np.nan]), ('source_sha256','stale')])
def test_reproduction_rejects_corrupted_physics_or_sources(key,value):
    reference = {'temperature':[310.], 'liquid_fraction':[.5], 'position_m':.001,
                 'stored_J':10., 'source_sha256':'original'}
    with pytest.raises(RuntimeError, match='reproduction'):
        verify_reproduction({**reference,key:value}, reference)
