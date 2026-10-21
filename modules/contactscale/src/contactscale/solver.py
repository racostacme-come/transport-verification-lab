"""Conservative 1-D finite-volume heat flow with insulated outer boundaries."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.linalg import eigh_tridiagonal


@dataclass(frozen=True)
class Slab:
    """Length [m], conductivity [W/(m K)], volumetric heat capacity [J/(m³ K)]."""

    length: float
    conductivity: float
    capacity: float

    def __post_init__(self):
        if not all(
            np.isfinite(x) and x > 0 for x in (self.length, self.conductivity, self.capacity)
        ):
            raise ValueError("Slab properties must be finite and strictly positive")

    @property
    def effusivity(self):
        return np.sqrt(self.conductivity * self.capacity)


class ThermalChain:
    """Assemble a thermal chain and cache its symmetric eigendecomposition.

    All arrays are copied. A face contact resistance [m² K/W] is added to
    the two half-cell resistances. At least two cells are required.
    This dense modal approach uses O(N²) storage; intended for modest 1-D grids.
    """

    def __init__(self, widths, conductivity, capacity, contacts=None):
        arrays = [np.array(v, dtype=float, copy=True) for v in (widths, conductivity, capacity)]
        if any(a.ndim != 1 for a in arrays) or len(arrays[0]) < 2:
            raise ValueError("Property arrays must be one-dimensional with at least two cells")
        if any(a.shape != arrays[0].shape for a in arrays):
            raise ValueError("Property arrays must have matching shapes")
        if any(not np.all(np.isfinite(a)) or np.any(a <= 0) for a in arrays):
            raise ValueError("Cell properties must be finite and positive")
        self.widths, self.conductivity, self.capacity = arrays
        n = len(self.widths)
        contacts = (
            np.zeros(n - 1) if contacts is None else np.array(contacts, dtype=float, copy=True)
        )
        if contacts.shape != (n - 1,) or not np.all(np.isfinite(contacts)) or np.any(contacts < 0):
            raise ValueError("Contacts must contain N-1 finite nonnegative resistances")
        self.contacts = contacts
        self.mass = self.widths * self.capacity
        half = self.widths / (2 * self.conductivity)
        self.conductance = 1 / (half[:-1] + contacts + half[1:])
        self.diagonal = np.r_[self.conductance, 0] + np.r_[0, self.conductance]
        roots = np.sqrt(self.mass)
        self.rates, self.modes = eigh_tridiagonal(
            self.diagonal / self.mass,
            -self.conductance / (roots[:-1] * roots[1:]),
            lapack_driver="stev",
        )
        # The connected insulated chain has one analytically zero eigenvalue.
        if abs(self.rates[0]) > 1e-10 * max(1, self.rates[-1]):
            raise ArithmeticError("Unexpected nonzero constant mode")
        self.rates[0] = 0
        for value in vars(self).values():
            if isinstance(value, np.ndarray):
                value.flags.writeable = False

    def _temperature(self, value):
        value = np.asarray(value, dtype=float)
        if value.shape != self.mass.shape or not np.all(np.isfinite(value)):
            raise ValueError("Temperature must be a finite vector of length N")
        return value

    def evolve(self, initial, times):
        """Return temperatures shaped (len(times), N); times need not be sorted.

        Computes exp(-M^-1 K t) without time stepping. Roundoff and spatial
        discretization remain. Exactly returns supplied initial data at t=0.
        """
        initial = self._temperature(initial)
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not np.all(np.isfinite(times)) or np.any(times < 0):
            raise ValueError("Times must be a finite nonnegative one-dimensional array")
        mean = np.dot(self.mass, initial) / np.sum(self.mass)
        roots = np.sqrt(self.mass)
        coefficients = self.modes.T @ (roots * (initial - mean))
        coefficients[0] = 0  # The mean is propagated analytically, without exponential drift.
        result = (
            mean + ((np.exp(-times[:, None] * self.rates) * coefficients) @ self.modes.T) / roots
        )
        result[times == 0] = initial
        return result

    def flux(self, temperature):
        """Return left-to-right face fluxes [W/m²] for a single temperature vector."""
        temperature = self._temperature(temperature)
        return self.conductance * (temperature[:-1] - temperature[1:])

    def surface_temperatures(self, temperature, face):
        """Return left/right surface temperatures reconstructed at an internal face."""
        if (
            isinstance(face, bool)
            or not isinstance(face, Integral)
            or not 0 <= face < len(self.contacts)
        ):
            raise ValueError("Face must be an integer in [0, N-2]")
        temperature = self._temperature(temperature)
        q = self.flux(temperature)[face]
        left = temperature[face] - q * self.widths[face] / (2 * self.conductivity[face])
        right = temperature[face + 1] + q * self.widths[face + 1] / (
            2 * self.conductivity[face + 1]
        )
        return left, right


def two_slabs(left, right, cells, resistance):
    """Build equal cell counts per slab; return (model, interface face index)."""
    if isinstance(cells, bool) or not isinstance(cells, Integral) or cells < 1:
        raise ValueError("cells must be a positive integer")
    if not np.isfinite(resistance) or resistance < 0:
        raise ValueError("Resistance must be finite and nonnegative")
    widths = np.repeat([left.length / cells, right.length / cells], cells)
    conductivity = np.repeat([left.conductivity, right.conductivity], cells)
    capacity = np.repeat([left.capacity, right.capacity], cells)
    contacts = np.zeros(2 * cells - 1)
    contacts[cells - 1] = resistance
    return ThermalChain(widths, conductivity, capacity, contacts), cells - 1
