import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad
from scipy.linalg import expm

from contactscale import Slab, ThermalChain, contact_flux, contact_heat, contact_scale, two_slabs


@pytest.mark.parametrize("time", [0, 1e-14, 1e-6, 0.01, 0.1024, 1, 100, 1e8])
def test_heat_against_independent_quadrature(time):
    # Substitute t*u² to remove the sqrt(t) derivative singularity at the origin.
    expected = quad(
        lambda u: 2 * time * u * contact_flux(time * u * u, 200, 8000, 2000, 2e-4),
        0,
        1,
        epsabs=1e-9,
        epsrel=1e-11,
    )[0]
    assert_allclose(contact_heat(time, 200, 8000, 2000, 2e-4), expected, rtol=2e-10, atol=1e-12)


def test_limits_and_crossover():
    assert contact_scale(8000, 2000, 2e-4) == pytest.approx(0.1024)
    assert contact_flux(0, 200, 8000, 2000, 2e-4) == 1e6
    assert contact_heat(0, 200, 8000, 2000, 0) == 0
    assert_allclose(
        contact_flux(1, 200, 8000, 2000, 1e-10), contact_flux(1, 200, 8000, 2000, 0), rtol=1e-12
    )
    assert_allclose(
        contact_heat(1, 200, 8000, 2000, 1e-10), contact_heat(1, 200, 8000, 2000, 0), rtol=1e-6
    )
    assert_allclose(
        contact_flux([0, 1], -200, 8000, 2000, 2e-4), -contact_flux([0, 1], 200, 8000, 2000, 2e-4)
    )


@pytest.mark.parametrize("resistance", [0, 0.2, 10])
def test_two_cell_exact(resistance):
    model = ThermalChain([0.2, 0.3], [2, 4], [3, 5], [resistance])
    t = np.array([0, 0.1, 1, 10])
    initial = np.array([7, 2])
    mass = model.mass
    mean = mass @ initial / mass.sum()
    decay = np.exp(-model.conductance[0] * (1 / mass).sum() * t)
    expected = mean + decay[:, None] * (initial - mean)
    assert_allclose(model.evolve(initial, t), expected, atol=2e-14)


def test_dense_exponential_energy_maximum_and_jump():
    model = ThermalChain([0.1, 0.2, 0.3, 0.1], [1, 4, 2, 3], [3, 2, 7, 2], [0, 0.3, 0.1])
    initial = np.array([700, 300, 550, 310])
    k = np.diag(model.diagonal) - np.diag(model.conductance, 1) - np.diag(model.conductance, -1)
    times = [0, 0.001, 0.2, 2, 100]
    result = model.evolve(initial, times)
    for t, temp in zip(times, result, strict=True):
        assert_allclose(temp, expm(-k / model.mass[:, None] * t) @ initial, atol=1e-9)
        assert_allclose(model.mass @ temp, model.mass @ initial, atol=1e-10)
        assert temp.min() >= initial.min() - 1e-10
        assert temp.max() <= initial.max() + 1e-10
        left, right = model.surface_temperatures(temp, 1)
        assert_allclose(left - right, model.flux(temp)[1] * model.contacts[1], atol=1e-12)
        assert_allclose(
            -temp @ k @ temp, -np.sum(model.conductance * np.diff(temp) ** 2), atol=1e-8
        )
    assert np.all(np.diff((result**2) @ model.mass) <= 1e-8)
    assert_allclose(model.evolve(np.ones(4) * 300, times), 300, atol=0)


def test_smooth_continuum_refinement():
    errors = []
    for n in (10, 20, 40):
        h = 1 / n
        model = ThermalChain(np.full(n, h), np.ones(n), np.ones(n))
        # Exact cell averages of a Neumann cosine, not point values.
        edges = np.linspace(0, 1, n + 1)
        initial = np.diff(np.sin(np.pi * edges)) / (np.pi * h)
        exact = initial * np.exp(-(np.pi**2) * 0.1)
        errors.append(np.linalg.norm(model.evolve(initial, [0.1])[0] - exact) / np.sqrt(n))
    assert np.all(np.log2(np.array(errors[:-1]) / errors[1:]) > 1.98)


@pytest.mark.parametrize("bad", [-1, np.inf, np.nan])
def test_bad_analytical_inputs(bad):
    with pytest.raises(ValueError):
        contact_flux(bad, 200, 8000, 2000, 1e-3)
    with pytest.raises(ValueError):
        contact_heat(1, 200, 8000, 2000, bad)
    with pytest.raises(ValueError):
        Slab(bad, 1, 1)


def test_bad_shapes_and_times():
    with pytest.raises(ValueError):
        ThermalChain([1, 2], [1], [1, 2])
    with pytest.raises(ValueError):
        ThermalChain([1, 2], [1, 1], [1, 2], [-1])
    model = ThermalChain([1, 2], [1, 1], [1, 2])
    with pytest.raises(ValueError):
        model.evolve([1, 2], -1)
    with pytest.raises(ValueError):
        model.evolve([1, np.nan], [1])
    with pytest.raises(ValueError):
        model.surface_temperatures([1, 2], 1)
    with pytest.raises(ValueError):
        contact_flux(0, 200, 8000, 2000, 0)
    for cells in (0, 1.5, True):
        with pytest.raises(ValueError):
            two_slabs(Slab(1, 1, 1), Slab(1, 1, 1), cells, 0)
