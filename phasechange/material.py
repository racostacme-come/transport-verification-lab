"""Monotone piecewise-linear volumetric enthalpy with an explicit inverse."""
from dataclasses import dataclass, asdict
import math
import numpy as np


def finite_array(value, name):
    result = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def finite_result(value, name):
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} exceeds floating-point range")
    return value


@dataclass(frozen=True)
class Material:
    """Constant rho [kg/m3], c [J/kg/K], k [W/m/K], L [J/kg], T [K].

    H= rho*c*(T-reference) + rho*L*f(T) in J/m3. The phase fraction is
    constitutively saturated outside [solidus, liquidus], not a state repair.
    A zero latent heat is supported; an isothermal phase jump is not.
    """
    density: float
    heat_capacity: float
    conductivity: float
    latent_heat: float
    solidus: float
    liquidus: float
    reference: float = 0.0

    def __post_init__(self):
        for name, value in asdict(self).items():
            if not math.isfinite(value):
                raise ValueError(f"{name} must be finite")
        if min(self.density, self.heat_capacity, self.conductivity) <= 0:
            raise ValueError("density, heat capacity and conductivity must be positive")
        if self.latent_heat < 0 or self.solidus < 0 or self.reference < 0:
            raise ValueError("latent heat and absolute temperatures must be nonnegative")
        if self.liquidus <= self.solidus:
            raise ValueError("liquidus must exceed solidus")
        derived = (self.capacity, self.latent_density, self.mushy_capacity,
                   self.capacity*(self.solidus-self.reference),
                   self.capacity*(self.liquidus-self.reference)+self.latent_density)
        if not all(math.isfinite(x) for x in derived) or self.capacity <= 0:
            raise ValueError("derived enthalpy scales must be finite and capacity positive")

    @property
    def capacity(self):
        return self.density * self.heat_capacity

    @property
    def latent_density(self):
        return self.density * self.latent_heat

    @property
    def mushy_capacity(self):
        return self.capacity + self.latent_density / (self.liquidus-self.solidus)

    def _temperature(self, temperature):
        value = finite_array(temperature, "temperature")
        if np.any(value < 0):
            raise ValueError("temperature must be nonnegative kelvin")
        return value

    def liquid_fraction(self, temperature):
        value = self._temperature(temperature)
        with np.errstate(over="ignore"):
            return np.clip((value-self.solidus)/(self.liquidus-self.solidus), 0., 1.)

    def enthalpy(self, temperature):
        value = self._temperature(temperature)
        with np.errstate(over="ignore", invalid="ignore"):
            result = self.capacity*(value-self.reference) + self.latent_density*self.liquid_fraction(value)
        return finite_result(result, "enthalpy")

    def temperature(self, enthalpy):
        """Exact branch inverse; values implying a negative Kelvin state fail."""
        value = finite_array(enthalpy, "enthalpy")
        hs = self.capacity*(self.solidus-self.reference)
        hl = self.capacity*(self.liquidus-self.reference)+self.latent_density
        # Endpoint-relative form limits cancellation in the phase interval.
        with np.errstate(over="ignore", invalid="ignore"):
            result = np.where(value <= hs, self.reference+value/self.capacity,
                              np.where(value < hl, self.solidus+(value-hs)/self.mushy_capacity,
                                       self.liquidus+(value-hl)/self.capacity))
        finite_result(result, "temperature")
        return self._temperature(result)

    def tangent(self, temperature, *, side="right"):
        """dH/dT [J/m3/K]; explicit one-sided value at either phase endpoint."""
        value = self._temperature(temperature)
        if side not in ("left", "right"):
            raise ValueError("side must be left or right")
        inside = ((value >= self.solidus) & (value < self.liquidus) if side == "right"
                  else (value > self.solidus) & (value <= self.liquidus))
        return np.where(inside, self.mushy_capacity, self.capacity)
