import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import quad

from contactscale import contact_flux, contact_heat


@pytest.mark.parametrize("time", [1e-6, 0.001, 0.1, 1, 100])
def test_continuum_surface_impedance(time):
    """The half-space convolution plus contact jump must recover the initial gap."""
    args = (200, 8000, 2000, 2e-4)
    b = 1 / 8000 + 1 / 2000
    integral = (
        2
        * np.sqrt(time)
        * quad(
            lambda v: contact_flux(time * (1 - v * v), *args),
            0,
            1,
            epsabs=1e-7,
            epsrel=1e-10,
        )[0]
    )
    assert_allclose(
        2e-4 * contact_flux(time, *args) + b / np.sqrt(np.pi) * integral, 200, rtol=1e-10
    )


@pytest.mark.parametrize("z", [0.099999, 0.1, 0.100001, 1, 10])
def test_integrated_heat_derivative_and_series_join(z):
    time = (z * 2e-4 / 0.000625) ** 2
    dt = time * 1e-4
    args = (200, 8000, 2000, 2e-4)
    derivative = (contact_heat(time + dt, *args) - contact_heat(time - dt, *args)) / (2 * dt)
    assert_allclose(derivative, contact_flux(time, *args), rtol=2e-8)
