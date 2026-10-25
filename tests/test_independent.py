"""Independent two-cell and Fourier references; no shared eigensolver oracle."""
import numpy as np
import pytest
from contactscale.solver import ThermalChain
from heatstart.solver import HeatGrid, evolve, exact_discrete
from verification import experiment


@pytest.mark.parametrize('method', ['modal', 'be', 'cn', 'rannacher'])
def test_two_unequal_cells_have_analytic_exchange(method):
    # M=(0.6, 2.8), g=1/(0.3/4 + 0.7/10); unequal volumes/materials.
    grid = HeatGrid([0, 0.3, 1], [2, 5], [2, 4])
    chain = ThermalChain([0.3, 0.7], [2, 5], [2, 4])
    mass = np.array([0.6, 2.8])
    conductance = 1 / (0.3 / 4 + 0.7 / 10)
    rate = conductance * (1 / mass[0] + 1 / mass[1])
    initial = np.array([3., -1.])
    mean = mass @ initial / mass.sum()
    dt, steps = 0.01, 7
    if method == 'modal':
        factor = np.exp(-rate * dt * steps)
        expected = mean + factor * (initial - mean)
        np.testing.assert_allclose(chain.evolve(initial, [dt * steps])[0], expected, atol=2e-13)
        actual = exact_discrete(grid, initial, dt * steps)
    else:
        be = 1 / (1 + rate * dt)
        cn = (1 - rate * dt / 2) / (1 + rate * dt / 2)
        factor = {'be': be**steps, 'cn': cn**steps,
                  'rannacher': (1 + rate * dt / 2)**(-4) * cn**(steps - 2)}[method]
        expected = mean + factor * (initial - mean)
        actual = evolve(grid, initial, dt, steps, method=method)[-1]
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)


def test_smooth_mode_separates_temporal_and_continuum_orders():
    result = experiment()
    for row in result['continuum']:
        assert row['heatstart_analytic_max_error'] < 1e-11
        assert row['contactscale_analytic_max_error'] < 1e-11
        if 'observed_order' in row:
            assert 1.98 < row['observed_order'] < 2.02
    for row in result['temporal']:
        if row['steps'] == 64:
            target = 1 if row['method'] == 'be' else 2
            assert abs(row['observed_order'] - target) < 0.04


def test_local_explicit_positivity_bound_is_sharp():
    grid = HeatGrid([0, .08, .3, .65, 1], [1, 3, .5, 2], [2, 1, 3, 1])
    limiting = int(np.argmin(grid.mass / grid.diagonal))
    initial = np.eye(4)[limiting]
    for fraction in (0.99, 1.01):
        result = initial - fraction * grid.explicit_bound * grid.stiffness_action(initial) / grid.mass
        assert abs(grid.mass @ (result - initial)) < 1e-15
        if fraction < 1:
            assert result.min() >= 0
        else:
            assert result[limiting] == pytest.approx(-.01)


@pytest.mark.parametrize('method', ['be', 'cn', 'rannacher'])
def test_nonuniform_grid_preserves_heat_and_dissipates_quadratic_energy(method):
    grid = HeatGrid([0, .08, .3, .65, 1], [1, 3, .5, 2], [2, 1, 3, 1])
    initial = np.array([1., -2., 4., 0.])
    history = evolve(grid, initial, .1, 12, method=method)
    np.testing.assert_allclose(history @ grid.mass, grid.mass @ initial, rtol=1e-12, atol=1e-12)
    assert np.max(np.diff(.5 * history**2 @ grid.mass)) < 1e-13
    # A constant offset must remain unchanged by the evolution.
    shifted = evolve(grid, initial + 7, .1, 12, method=method)
    np.testing.assert_allclose(shifted - history, 7, rtol=1e-12, atol=1e-12)
