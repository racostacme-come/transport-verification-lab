"""Reproduce conservative slab checks: python -m phasechange.slab_checks."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

import numpy as np
from .finite_volume import Slab, advance
from .material import Material

ROOT = Path(__file__).resolve().parents[1]


def experiment():
    params = json.loads((ROOT/'phasechange/benchmarks.json').read_text())['material']
    material = Material(**params)
    slab = Slab(length=0.012, cells=24, area=0.003)
    initial = np.full(slab.cells, 305.0)
    temperature = initial
    rows = []
    for wall, power in ((335., 2.7e6), (295., -2.7e6)):
        for _ in range(5):
            step = advance(material, slab, temperature, 10., left=wall, source=power)
            rows.append(dict(wall_K=wall, source_W_m3=power, **asdict(step)))
            temperature = step.temperature

    # A cosine is an eigenvector of the insulated cell-centred diffusion
    # operator. Backward Euler has an exact algebraic amplification factor.
    linear = replace(material, latent_heat=0.)
    coordinate = (np.arange(slab.cells)+0.5)/slab.cells
    mode = 8*np.cos(np.pi*coordinate)
    dt = 12.
    rate = 4*linear.conductivity/linear.capacity*(np.sin(np.pi/(2*slab.cells))/slab.dx)**2
    linear_step = advance(linear, slab, 310+mode, dt)
    mode_error = float(np.max(np.abs(linear_step.temperature-(310+mode/(1+dt*rate)))))
    result = {
        'slab': asdict(slab), 'initial_K': initial.tolist(), 'dt_s': 10., 'steps': rows,
        'insulated_cosine': {'dt_s': dt, 'error_K': mode_error, 'decay_rate_per_s': rate},
        'max_scaled_residual': max(step['scaled_residual'] for step in rows),
        'max_scaled_energy_defect': max(step['ledger']['scaled_defect'] for step in rows),
    }
    files = ('finite_volume.py', 'material.py', 'slab_checks.py', 'benchmarks.json', '__init__.py')
    result['source_sha256'] = {f'phasechange/{name}': hashlib.sha256(
        (ROOT/'phasechange'/name).read_bytes()).hexdigest() for name in files}
    if not (mode_error <= 1e-10 and result['max_scaled_residual'] <= 1e-10
            and result['max_scaled_energy_defect'] <= 1e-9):
        raise RuntimeError("slab checks failed")
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results/enthalpy-slab.json')
    args = parser.parse_args()
    data = experiment()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print({k: data[k] for k in ('insulated_cosine', 'max_scaled_residual', 'max_scaled_energy_defect')})
