"""Conservative cell-centred backward Euler for a uniform 1-D slab (SI units)."""
from dataclasses import dataclass
from math import fsum, isfinite
from numbers import Integral

import numpy as np
from scipy.linalg import solveh_banded
from scipy.optimize import brentq

from .material import Material, finite_array, finite_result


@dataclass(frozen=True)
class Slab:
    length: float
    cells: int
    area: float = 1.0

    def __post_init__(self):
        if isinstance(self.cells, bool) or not isinstance(self.cells, Integral) or self.cells < 1:
            raise ValueError("cells must be a positive integer")
        if not all(isfinite(x) and x > 0 for x in (self.length, self.area)):
            raise ValueError("length and area must be finite and positive")

    @property
    def dx(self):
        return self.length / self.cells


@dataclass(frozen=True)
class Ledger:
    stored_J: float
    left_J: float
    right_J: float
    source_J: float
    defect_J: float
    scaled_defect: float


@dataclass(frozen=True)
class Step:
    temperature: tuple
    enthalpy: tuple
    liquid_fraction: tuple
    iterations: int
    scaled_residual: float
    ledger: Ledger


def _field(values, cells, name):
    field = finite_array(values, name)
    if field.ndim == 0:
        return np.full(cells, field.item())
    if field.shape != (cells,):
        raise ValueError(f"{name} must be scalar or contain one value per cell")
    return field


def _ledger(mat, slab, old, new, dt, source, left, right, scale):
    # Integrate the capacity over the traversed temperature interval. This
    # reconstruction does not use the residual, its face array or H(new)-H(old).
    low = np.maximum(np.minimum(old, new), mat.solidus)
    high = np.minimum(np.maximum(old, new), mat.liquidus)
    span = np.sign(new-old) * np.maximum(high-low, 0.0)
    volume = slab.area * slab.dx
    stored = volume * fsum(mat.capacity*(new-old) +
                           mat.latent_density*(span/(mat.liquidus-mat.solidus)))
    def wall_heat(wall, surface_cell):
        if wall is None:
            return 0.0
        return dt * slab.area * mat.conductivity * (wall-surface_cell) / (0.5*slab.dx)

    lhs = wall_heat(left, new[0])
    rhs = wall_heat(right, new[-1])
    added = dt * volume * fsum(source)
    defect = fsum((stored, -lhs, -rhs, -added))
    return Ledger(stored, lhs, rhs, added, defect, abs(defect)/(slab.cells*scale))


def advance(material: Material, slab: Slab, previous, dt, *, source=0.0,
            left=None, right=None, max_iterations=60):
    """Advance one step; walls are prescribed kelvin or None for insulation.

    previous [K] and source [W/m3] can be scalar or cell arrays. Positive
    ledger entries are heat entering the slab. Failure returns no partial step.
    """
    if not isfinite(dt) or dt <= 0:
        raise ValueError("dt must be positive and finite")
    if not isinstance(max_iterations, Integral) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if any(w is not None and (not isfinite(w) or w < 0) for w in (left, right)):
        raise ValueError("walls must be finite nonnegative kelvin or None")

    n = slab.cells
    old = _field(previous, n, "previous")
    source = _field(source, n, "source")
    old_fraction = material.liquid_fraction(old)
    width = material.liquidus - material.solidus
    volume = slab.area * slab.dx
    face = dt * material.conductivity * slab.area / slab.dx
    scale = volume * (material.capacity * width + material.latent_density)
    if min(face, scale) <= 0 or not all(isfinite(x) for x in (face, scale, n*scale)):
        raise ValueError("unrepresentable face or storage scale")

    def residual(temp):
        # Extend the solid branch through zero only for Newton trials.
        # Work with increments to avoid cancellation against an arbitrary Tref.
        fraction = np.clip((temp-material.solidus)/width, 0.0, 1.0)
        change = material.capacity*(temp-old) + material.latent_density*(fraction-old_fraction)
        values = volume*(change-dt*source)
        exchange = face*np.diff(temp)
        values[:-1] -= exchange
        values[1:] += exchange
        if left is not None:
            values[0] += 2*face*(temp[0]-left)
        if right is not None:
            values[-1] += 2*face*(temp[-1]-right)
        return finite_result(values, "energy residual")

    diagonal = np.zeros(n)
    diagonal[:-1] += face
    diagonal[1:] += face
    if left is not None:
        diagonal[0] += 2*face
    if right is not None:
        diagonal[-1] += 2*face
    bands = np.zeros((2, n))
    bands[1, :-1] = -face
    temp = old.copy()
    for iteration in range(max_iterations+1):
        values = residual(temp)
        error = float(np.max(np.abs(values))/scale)
        if error <= 1e-10:
            enthalpy = material.enthalpy(temp)  # validates accepted kelvin state
            ledger = _ledger(material, slab, old, temp, dt, source, left, right, scale)
            if not isfinite(ledger.scaled_defect) or ledger.scaled_defect > 1e-9:
                raise RuntimeError("global energy ledger failed")
            return Step(tuple(temp), tuple(enthalpy), tuple(material.liquid_fraction(temp)),
                        iteration, error, ledger)
        if iteration == max_iterations:
            break
        mushy = (temp >= material.solidus) & (temp < material.liquidus)
        capacity = np.where(mushy, material.mushy_capacity, material.capacity)
        bands[0] = volume*capacity + diagonal
        direction = (solveh_banded(bands, -values, lower=True) if n > 1
                     else -values/bands[0])
        # R is the gradient of a strictly convex potential. Its directional
        # derivative increases monotonically along the Newton search line.
        def slope(alpha):
            return float(np.dot(direction, residual(temp + alpha*direction)))

        alpha = 1.0
        if slope(1.0) > 0:
            alpha = brentq(slope, 0.0, 1.0, xtol=1e-15)
        temp = temp + alpha*direction
    raise RuntimeError(f"iteration limit {max_iterations} reached; scaled residual {error:.3e}")
