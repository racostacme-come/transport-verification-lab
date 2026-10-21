"""Semi-infinite contact: stable evaluation without exp(z*z) overflow."""

import numpy as np
from scipy.special import erfcx


def _parameters(time, delta, effusivity_left, effusivity_right, resistance):
    values = np.asarray([delta, effusivity_left, effusivity_right, resistance], dtype=float)
    if not np.all(np.isfinite(values)) or min(values[1:3]) <= 0 or resistance < 0:
        raise ValueError("Finite parameters, positive effusivities and nonnegative R required")
    time = np.asarray(time, dtype=float)
    if not np.all(np.isfinite(time)) or np.any(time < 0):
        raise ValueError("Time must be finite and nonnegative")
    return time, 1.0 / effusivity_left + 1.0 / effusivity_right


def contact_scale(effusivity_left, effusivity_right, resistance):
    """Return crossover time (R / (1/e_left + 1/e_right))**2, in seconds."""
    _, b = _parameters(0, 0, effusivity_left, effusivity_right, resistance)
    return (resistance / b) ** 2


def contact_flux(time, delta, effusivity_left, effusivity_right, resistance):
    """Flux [W/m²], positive from initially hotter left body to right.

    For R=0, time must be strictly positive: the initial flux is singular.
    Scalars and arbitrary-shaped time arrays are supported.
    """
    time, b = _parameters(time, delta, effusivity_left, effusivity_right, resistance)
    if resistance == 0:
        if np.any(time == 0):
            raise ValueError("Perfect-contact flux is undefined at time zero")
        return delta / (b * np.sqrt(np.pi * time))
    return (delta / resistance) * erfcx(b * np.sqrt(time) / resistance)


def contact_heat(time, delta, effusivity_left, effusivity_right, resistance):
    """Integrated contact heat [J/m²] from zero to time.

    A convergent Taylor series avoids cancellation for z=b*sqrt(t)/R<0.1.
    No numerical time integration is used.
    """
    time, b = _parameters(time, delta, effusivity_left, effusivity_right, resistance)
    if resistance == 0:
        return 2 * delta * np.sqrt(time / np.pi) / b
    z = b * np.sqrt(time) / resistance
    small = z < 0.1
    # Integral of erfcx(z*sqrt(u)) on 0<=u<=1, by its power series.
    coefficients = [1.0, -2 / np.sqrt(np.pi)]
    for n in range(2, 20):
        coefficients.append(2 * coefficients[n - 2] / n)
    series = np.polynomial.polynomial.polyval(
        np.where(small, z, 0), [2 * c / (n + 2) for n, c in enumerate(coefficients)]
    )
    safe_z = np.where(small, 1, z)
    large = delta / b * np.sqrt(time) * (2 / np.sqrt(np.pi) + (erfcx(safe_z) - 1) / safe_z)
    return np.where(small, delta / resistance * time * series, large)
