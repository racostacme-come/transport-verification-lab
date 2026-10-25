"""Cross-check independent thermal operators and isolate time integration error."""
import numpy as np
from contactscale.solver import ThermalChain
from heatstart.solver import HeatGrid, evolve, exact_discrete

def experiment():
    n = 40
    edges = np.linspace(0, 1, n+1)
    widths = np.diff(edges)
    k = np.where(edges[:-1] < 0.5, 1.0, 3.0)
    capacity = np.where(edges[:-1] < 0.5, 2.0, 1.0)
    grid = HeatGrid(edges, k, capacity)
    chain = ThermalChain(widths, k, capacity)
    initial = np.where(grid.centres < 0.5, 1.0, 0.0)
    final_time = 0.02
    reference = chain.evolve(initial, [final_time])[0]
    independent = exact_discrete(grid, initial, final_time)
    rows = []
    for steps in (4, 8, 16, 32):
        for method in ("be", "cn", "rannacher"):
            result = evolve(grid, initial, final_time/steps, steps, method=method)[-1]
            rows.append(dict(steps=steps, method=method,
                mass_weighted_error=float(np.sqrt(np.dot(grid.mass, (result-reference)**2))),
                heat_balance_error=float(abs(grid.mass @ (result-initial))),
                minimum_temperature=float(result.min())))
    return {"experiment": "heterogeneous-insulated-heat-chain", "rows": rows,
        "mass_operator_max_difference": float(np.max(abs(grid.mass-chain.mass))),
        "conductance_max_difference": float(np.max(abs(grid.conductance-chain.conductance))),
        "modal_reference_max_difference": float(np.max(abs(reference-independent))),
        "scope": "Zero contact resistance is the shared model. Error is against a fixed discrete grid, not a continuum error."}
