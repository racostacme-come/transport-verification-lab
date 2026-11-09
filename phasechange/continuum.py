"""Independent analytical oracles; no constitutive or solver imports.

One-phase melting with constant properties, a solid initially at Tm and an
isothermal hot wall. The solid remains at Tm. This sharp-interface reference
is distinct from the numerical finite melting interval.
"""
import math
import numpy as np
from scipy.special import erf


def similarity_parameter(stefan_number):
    """Solve lambda*erf(lambda)=Ste*exp(-lambda**2)/sqrt(pi)."""
    if not math.isfinite(stefan_number) or stefan_number <= 0:
        raise ValueError('Stefan number must be finite and positive')
    low, high = 0., min(1., math.sqrt(stefan_number))
    target = math.log(stefan_number)-.5*math.log(math.pi)
    def log_product(value):
        return math.log(value)+math.log(math.erf(value))+value*value
    while log_product(high) < target:
        high *= 2
    # Independent bisection; the left side increases and the right decreases.
    for _ in range(110):
        mid = (low+high)/2
        if log_product(mid) < target:
            low = mid
        else:
            high = mid
    return (low+high)/2


def stefan_solution(x, time, *, diffusivity, stefan_number, melting, wall):
    """Return sharp front [m] and point temperatures [K] at t>0."""
    x = np.asarray(x, dtype=float)
    if (not all(math.isfinite(v) for v in (time, diffusivity, melting, wall))
            or time <= 0 or diffusivity <= 0 or melting < 0 or wall <= melting
            or not np.all(np.isfinite(x)) or np.any(x < 0)):
        raise ValueError('invalid one-phase melting inputs')
    lam = similarity_parameter(stefan_number)
    scale = 2*math.sqrt(diffusivity*time)
    front = lam*scale
    temperature = np.where(x < front,
                           wall-(wall-melting)*erf(x/scale)/math.erf(lam), melting)
    return front, temperature


def cosine_reference(cells, length, diffusivity, time, *, discrete=False):
    """Unit cosine cell averages: continuum or exact FV semidiscrete evolution.

    Insulated ends. Initial data are exact cell averages in both models.
    discrete=True changes only the Laplacian eigenvalue, not time integration.
    """
    if (isinstance(cells, bool) or not isinstance(cells, (int, np.integer)) or cells < 2
            or not all(math.isfinite(v) for v in (length, diffusivity, time))
            or min(length, diffusivity) <= 0 or time < 0):
        raise ValueError('invalid cosine reference inputs')
    dx = length/cells
    rate = (4*diffusivity*(math.sin(math.pi/(2*cells))/dx)**2 if discrete
            else diffusivity*(math.pi/length)**2)
    averages = np.sinc(1/(2*cells))*np.cos(math.pi*(np.arange(cells)+.5)/cells)
    return averages*math.exp(-rate*time)
