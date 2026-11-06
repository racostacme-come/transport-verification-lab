"""Independent discrete oracles, phase transitions, and solver failure paths."""
from dataclasses import replace
import json
import numpy as np
import pytest
from phasechange.finite_volume import Slab, advance
from phasechange.material import Material
from phasechange.reference import event_temperature
from phasechange.slab_checks import ROOT, experiment

MAT = Material(2., 3., 1.5, 12., 10., 14.)
SLAB = Slab(2., 7, 0.3)


@pytest.mark.parametrize('initial,power,dt', [(8,12,10), (24,-12,10), (10,-12,1),
                                             (14,12,1), (14,-12,1), (10,12,1), (12,0,2)])
@pytest.mark.parametrize('cells', [1,7])
def test_uniform_crossings_and_endpoint_directions(initial, power, dt, cells):
    mesh = replace(SLAB, cells=cells)
    step = advance(MAT, mesh, initial, dt, source=power)
    expected = event_temperature(initial, dt, power, density=2, heat_capacity=3,
                                 latent_heat=12, solidus=10, liquidus=14)
    np.testing.assert_allclose(step.temperature, expected, atol=1e-12, rtol=0)
    assert step.ledger.stored_J == pytest.approx(mesh.area*mesh.length*dt*power, abs=1e-12)


@pytest.mark.parametrize('cells', [1,2,23])
def test_dirichlet_linear_profile_crossing_mushy_layer_is_steady(cells):
    mesh = replace(SLAB, cells=cells)
    initial = 7+12*(np.arange(cells)+0.5)/cells
    step = advance(MAT, mesh, initial, 2.5, left=7., right=19.)
    flux = 2.5*MAT.conductivity*mesh.area*12/mesh.length
    np.testing.assert_allclose(initial, step.temperature, atol=1e-12, rtol=0)
    assert step.ledger.left_J == pytest.approx(-flux)
    assert step.ledger.right_J == pytest.approx(flux)


def test_mixed_boundary_nonuniform_source_recovers_manufactured_discrete_state():
    old = np.array([18., 9., 13., 12., 8., 16., 14.])
    new = old[::-1]
    left, dt = 20., 3.
    # Direct power balance assembled cell-by-cell, including the half-cell wall.
    volume = SLAB.dx*SLAB.area
    def h(t):
        return 6*t + 6*min(4., max(0., t-10.))
    powers = []
    for i,t in enumerate(new):
        west = new[i-1] if i else 2*left-t
        east = new[i+1] if i+1 < len(new) else t
        powers.append((h(t)-h(old[i]))/dt - MAT.conductivity*(west+east-2*t)/SLAB.dx**2)
    step = advance(MAT, SLAB, old, dt, source=powers, left=left)
    np.testing.assert_allclose(step.temperature, new, atol=1e-11, rtol=0)
    assert step.ledger.source_J == pytest.approx(dt*volume*sum(powers))
    assert abs(step.ledger.defect_J) < 1e-10


def test_iteration_exhaustion_preserves_input():
    initial = np.linspace(8, 16, SLAB.cells)
    original = initial.copy()
    with pytest.raises(RuntimeError, match='iteration limit'):
        advance(MAT, SLAB, initial, 2, left=95, source=10, max_iterations=1)
    np.testing.assert_array_equal(initial, original)


@pytest.mark.parametrize('args', [dict(left=-1), dict(source=[1,2]), dict(right=np.nan),
                                  dict(max_iterations=0), dict(source=-1e5)])
def test_invalid_inputs_or_unphysical_solution_fail(args):
    with pytest.raises(ValueError):
        advance(MAT, SLAB, 12, 1, **args)


def test_versioned_slab_history_reproduces():
    actual = experiment()
    recorded = json.loads((ROOT/'results/enthalpy-slab.json').read_text())
    assert actual['source_sha256'] == recorded['source_sha256']
    assert len(actual['steps']) == len(recorded['steps']) == 10
    for step, reference in zip(actual['steps'], recorded['steps']):
        np.testing.assert_allclose(step['temperature'], reference['temperature'], atol=1e-9, rtol=0)
        for name in ('stored_J','left_J','right_J','source_J','defect_J'):
            assert step['ledger'][name] == pytest.approx(reference['ledger'][name], abs=1e-8)
    assert min(actual['steps'][4]['liquid_fraction']) == 1
    assert min(actual['steps'][-1]['liquid_fraction']) == 0
