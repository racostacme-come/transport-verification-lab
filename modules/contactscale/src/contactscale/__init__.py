"""Linear thermal contact models, with SI quantities per unit interface area."""

from .analytical import contact_flux, contact_heat, contact_scale
from .solver import Slab, ThermalChain, two_slabs

__all__ = ["Slab", "ThermalChain", "contact_flux", "contact_heat", "contact_scale", "two_slabs"]
