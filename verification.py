"""Analytic references separating implementation, temporal and continuum error.

Insulated unit interval, unit area, constant k=1.7 and capacity c=2.3.
The smooth initial field is 1 + cos(pi*x); arrays are exact cell averages.
No eigenvalues/eigenvectors from either implementation define these references.
"""
import numpy as np
from contactscale.solver import ThermalChain
from heatstart.solver import HeatGrid, evolve, exact_discrete


def cosine_case(cells):
    edges = np.linspace(0.0, 1.0, cells + 1)
    grid = HeatGrid(edges, np.full(cells, 1.7), np.full(cells, 2.3))
    chain = ThermalChain(np.diff(edges), grid.conductivity, grid.capacity)
    mode = np.sinc(1 / (2 * cells)) * np.cos(np.pi * grid.centres)
    rate = (1.7 / 2.3) * 4 * cells**2 * np.sin(np.pi / (2 * cells))**2
    return grid, chain, mode, rate


def norm(grid, value):
    return float(np.sqrt(grid.mass @ np.square(value)))


def experiment():
    time = 0.1
    continuum = []
    for cells in (10, 20, 40, 80):
        grid, chain, mode, rate = cosine_case(cells)
        initial = 1 + mode
        analytic_discrete = 1 + mode * np.exp(-rate * time)
        analytic_continuum = 1 + mode * np.exp(-(1.7 / 2.3) * np.pi**2 * time)
        heat = exact_discrete(grid, initial, time)
        contact = chain.evolve(initial, [time])[0]
        error = norm(grid, contact - analytic_continuum)
        row = dict(cells=cells, continuum_mass_error=error,
                   heatstart_analytic_max_error=float(np.max(abs(heat - analytic_discrete))),
                   contactscale_analytic_max_error=float(np.max(abs(contact - analytic_discrete))))
        if continuum:
            row['observed_order'] = float(np.log2(continuum[-1]['continuum_mass_error'] / error))
        continuum.append(row)
    grid, chain, mode, rate = cosine_case(40)
    reference = 1 + mode * np.exp(-rate * time)
    temporal = []
    for method in ('be', 'cn', 'rannacher'):
        previous = None
        for steps in (8, 16, 32, 64):
            result = evolve(grid, 1 + mode, time / steps, steps, method=method)[-1]
            error = norm(grid, result - reference)
            row = dict(method=method, steps=steps, discrete_mass_error=error)
            if previous is not None:
                row['observed_order'] = float(np.log2(previous / error))
            previous = error
            temporal.append(row)
    return dict(experiment='analytic-insulated-cosine', final_time=time,
                conductivity=1.7, capacity=2.3, continuum=continuum, temporal=temporal,
                scope='Cell-average continuum reference; discrete cosine eigenmode derived analytically. Synthetic verification only.')


def write_table(data, path):
    lines = [r'\begin{center}\small', r'\begin{tabular}{rrrr}\toprule',
             r'Cells & Continuum error & Order & Analytic discrete defect\\\midrule']
    for row in data['continuum']:
        order = f"{row['observed_order']:.3f}" if 'observed_order' in row else '--'
        defect = max(row['heatstart_analytic_max_error'], row['contactscale_analytic_max_error'])
        lines.append(f"{row['cells']} & {row['continuum_mass_error']:.3e} & {order} & {defect:.3e}" + r'\\')
    lines += [r'\bottomrule\end{tabular}', r'\par Smooth-mode temporal orders (32 to 64 steps):']
    lines.append(', '.join(f"{r['method']}: {r['observed_order']:.3f}" for r in data['temporal'] if r['steps'] == 64) + '.')
    lines.append(r'\end{center}')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8', newline='\n')
