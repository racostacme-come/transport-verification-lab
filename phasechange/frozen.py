"""Deliberately frozen apparent capacity, with the true enthalpy ledger retained."""
from math import isfinite
import numpy as np
from scipy.linalg import solveh_banded

from .finite_volume import Step, _field, _ledger
from .material import finite_result


def advance_frozen(material, slab, previous, dt, *, source=0., left=None,
                   right=None, backend=None):
    """Solve C(T_old) delta T = dt*(diffusion(T_new) + source).

    The right-sided tangent is fixed for the entire step, including cooling
    from a phase endpoint. This choice is part of the comparator definition.
    Step.scaled_residual checks this *linear* equation, not enthalpy balance.
    The independent true-enthalpy ledger is reported even when it is large.
    Native mode reuses the C++ Jacobian at T_old and its LDL solver explicitly.
    """
    if not isfinite(dt) or dt <= 0:
        raise ValueError('dt must be positive and finite')
    if any(w is not None and (not isfinite(w) or w < 0) for w in (left, right)):
        raise ValueError('walls must be finite nonnegative kelvin or None')
    n = slab.cells
    old = _field(previous, n, 'previous')
    power = _field(source, n, 'source')
    capacity = material.tangent(old, side='right')  # also validates kelvin
    volume = slab.area*slab.dx
    face = dt*material.conductivity*slab.area/slab.dx
    scale = volume*(material.capacity*(material.liquidus-material.solidus)
                    + material.latent_density)
    if min(face, scale) <= 0 or not all(isfinite(x) for x in (face, scale, n*scale)):
        raise ValueError('unrepresentable face or storage scale')

    # Work in increments, avoiding cancellation of a large absolute T offset.
    diagonal = volume*capacity.copy()
    rhs = dt*volume*power.copy()
    exchange = face*np.diff(old)
    rhs[:-1] += exchange
    rhs[1:] -= exchange
    diagonal[:-1] += face
    diagonal[1:] += face
    for i, wall in ((0, left), (n-1, right)):
        if wall is not None:
            rhs[i] += 2*face*(wall-old[i])
            diagonal[i] += 2*face
    finite_result(diagonal, 'frozen diagonal')
    finite_result(rhs, 'frozen right hand side')
    if backend is None:
        bands = np.zeros((2, n))
        bands[0] = diagonal
        bands[1, :-1] = -face
        increment = solveh_banded(bands, rhs, lower=True) if n > 1 else rhs/diagonal
    else:
        residual, native_diagonal, off = backend.assemble(
            material, slab, old, old, dt, source=power, left=left, right=right,
            side='right')
        increment = backend.solve(native_diagonal, off, -residual)
    temp = old + increment
    enthalpy = material.enthalpy(temp)  # reject, never clip, negative kelvin

    # Reconstruct the linear balance from the accepted rounded temperature.
    delta = temp-old
    residual = diagonal*delta-rhs
    residual[:-1] -= face*delta[1:]
    residual[1:] -= face*delta[:-1]
    error = float(np.max(np.abs(finite_result(residual, 'frozen residual')))/scale)
    if error > 1e-10:
        raise RuntimeError(f'frozen linear residual failed: {error:.3e}')
    ledger = _ledger(material, slab, old, temp, dt, power, left, right, scale)
    if not all(isfinite(x) for x in vars(ledger).values()):
        raise ValueError('nonfinite frozen energy ledger')
    return Step(tuple(temp), tuple(enthalpy), tuple(material.liquid_fraction(temp)),
                1, error, ledger)
