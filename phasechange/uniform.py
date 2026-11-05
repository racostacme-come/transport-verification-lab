"""Spatially uniform, insulated heating: dH/dt=q with piecewise constant q."""
from dataclasses import dataclass
import math
import numpy as np
from .material import Material, finite_array


@dataclass(frozen=True)
class UniformHistory:
    """Immutable samples including t=0; heat and energy defect in J/m3."""
    time: tuple[float, ...]
    temperature: tuple[float, ...]
    enthalpy: tuple[float, ...]
    liquid_fraction: tuple[float, ...]
    supplied_heat: tuple[float, ...]
    energy_defect: tuple[float, ...]


def integrate_uniform(material: Material, initial_temperature, durations, powers):
    """Integrate exact energy increments, never a frozen tangent approximation.

    Durations [s] and volumetric powers [W/m3] must be equally sized 1D
    arrays. Positive power heats; negative power cools. Sampling controls no
    discretization error for this prescribed piecewise constant input.
    """
    initial = finite_array(initial_temperature, "initial temperature")
    if initial.ndim != 0:
        raise ValueError("initial temperature must be scalar")
    dt = finite_array(durations, "durations")
    q = finite_array(powers, "powers")
    if dt.ndim != 1 or q.ndim != 1 or dt.shape != q.shape or dt.size == 0:
        raise ValueError("durations and powers must be nonempty matching 1D arrays")
    if np.any(dt <= 0):
        raise ValueError("durations must be positive")
    h0 = float(material.enthalpy(initial))
    h = [h0]; times = [0.]; temp = [float(initial)]; supplied = [0.]
    increments = []
    for step, power in zip(dt, q):
        increment = float(step)*float(power)
        if not math.isfinite(increment):
            raise ValueError("heat increment exceeds floating-point range")
        increments.append(increment)
        heat = math.fsum(increments)
        state = h[-1]+increment
        temperature = float(material.temperature(state))
        times.append(times[-1]+float(step))
        if not math.isfinite(times[-1]) or times[-1] <= times[-2]:
            raise ValueError("time grid is not representable as strictly increasing")
        temp.append(temperature); h.append(state); supplied.append(heat)
    # Reevaluate H from returned temperatures, independently of stored state.
    reconstructed = material.enthalpy(temp)
    defect = reconstructed-h0-np.asarray(supplied)
    return UniformHistory(tuple(times), tuple(temp), tuple(h),
                          tuple(map(float, material.liquid_fraction(temp))),
                          tuple(supplied), tuple(map(float, defect)))
