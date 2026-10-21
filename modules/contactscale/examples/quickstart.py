"""Run from any directory after installing contactscale."""

import numpy as np

from contactscale import Slab, contact_flux, two_slabs

left, right = Slab(0.06, 16, 4e6), Slab(0.03, 2, 2e6)
model, interface = two_slabs(left, right, cells=80, resistance=2e-4)
temperature = model.evolve(np.repeat([500.0, 300.0], 80), [1])[0]
print("Finite-volume flux:", model.flux(temperature)[interface], "W/m²")
print("Half-space reference:", contact_flux(1, 200, left.effusivity, right.effusivity, 2e-4))
