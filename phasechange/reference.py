"""Independent scalar event traversal for a uniform prescribed heat input.

Uses neither Material.enthalpy/temperature/tangent nor the energy integrator.
It spends time to reach phase boundaries, then changes the local dT/dt.
"""
import math


def event_temperature(initial, duration, power, *, density, heat_capacity,
                      latent_heat, solidus, liquidus):
    values = (initial, duration, power, density, heat_capacity, latent_heat, solidus, liquidus)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("reference inputs must be finite")
    if min(initial, duration, latent_heat, solidus) < 0 or min(density, heat_capacity) <= 0 or liquidus <= solidus:
        raise ValueError("invalid reference state or material")
    if duration == 0 or power == 0:
        return float(initial)
    temperature = float(initial); remaining = float(duration)
    capacity = density*heat_capacity
    mushy = density*(heat_capacity+latent_heat/(liquidus-solidus))
    # A monotone segment crosses at most two endpoints.
    for _ in range(3):
        if power > 0:
            boundary, slope = ((solidus, capacity) if temperature < solidus else
                               (liquidus, mushy) if temperature < liquidus else (math.inf, capacity))
        else:
            boundary, slope = ((liquidus, capacity) if temperature > liquidus else
                               (solidus, mushy) if temperature > solidus else (-math.inf, capacity))
        crossing_time = (boundary-temperature)*slope/power
        if remaining <= crossing_time:
            result = temperature + power*remaining/slope
            if not math.isfinite(result) or result < 0:
                raise ValueError("reference leaves representable Kelvin range")
            return result
        remaining -= crossing_time
        temperature = boundary
    raise RuntimeError("reference traversal exceeded two phase boundaries")
