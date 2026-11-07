"""Matched conservative/frozen experiments with source-bound, signed ledgers."""
from dataclasses import asdict, replace
from math import fsum, isfinite
import hashlib
import json
from pathlib import Path
import numpy as np

from .finite_volume import Slab, advance
from .frozen import advance_frozen
from .material import Material, finite_array
from .reference import event_temperature

ROOT = Path(__file__).resolve().parents[1]
METHODS = {'enthalpy': advance, 'frozen': advance_frozen}


def verify_reproduction(actual, recorded):
    """Fail on stale inputs or changed physical histories, on either backend.

    Temperatures: 2e-9 K absolute; fractions/front metres: 1e-9 absolute.
    Dimensional enthalpy/energy: 2e-7 absolute + 2e-9 relative (cancellation).
    Residual/closure maxima are gated during experiment(), not forced to have
    identical roundoff or iteration counts across linear algebra libraries.
    """
    def visit(a, b, key=''):
        if key in ('backend', 'iterations', 'scaled_residual',
                   'max_scaled_equation_residual', 'max_scaled_ledger_closure'):
            return
        if isinstance(b, dict):
            if not isinstance(a, dict) or a.keys() != b.keys():
                raise RuntimeError(f'reproduction structure mismatch: {key}')
            for k in b:
                visit(a[k], b[k], k)
        elif isinstance(b, (list, tuple)):
            if len(a) != len(b):
                raise RuntimeError(f'reproduction length mismatch: {key}')
            for x,y in zip(a,b):
                visit(x,y,key)
        elif isinstance(b, (float, int)):
            energy = key == 'enthalpy' or key.endswith('_J')
            atol, rtol = (2e-7, 2e-9) if energy else (1e-9, 0.)
            if key == 'temperature' or key.endswith('_K'):
                atol = 2e-9
            if not np.isfinite(a) or not np.isclose(a,b,atol=atol,rtol=rtol):
                raise RuntimeError(f'reproduction numerical mismatch: {key}: {a} != {b}')
        elif a != b:
            raise RuntimeError(f'reproduction metadata mismatch: {key}')
    visit(actual, recorded)


def front_position(fraction, slab):
    """Locate a unique decreasing f=0.5 crossing between cell centres.

    No wall extrapolation or silently chosen crossing. A front outside the
    centre interval, a plateau at 0.5, or a nonmonotone profile is unresolved.
    """
    f = finite_array(fraction, 'fraction')
    if f.shape != (slab.cells,) or np.any((f < 0) | (f > 1)):
        raise ValueError('fraction must contain one value in [0,1] per cell')
    if np.any(np.diff(f) > 0):
        return {'position_m': None, 'status': 'nonmonotone'}
    centers = (np.arange(slab.cells)+.5)*slab.dx
    exact = np.flatnonzero(f == .5)
    if len(exact) > 1:
        return {'position_m': None, 'status': 'plateau'}
    if len(exact) == 1:
        return {'position_m': float(centers[exact[0]]), 'status': 'resolved'}
    indices = np.flatnonzero((f[:-1] > .5) & (f[1:] < .5))
    if len(indices) == 0:
        return {'position_m': None, 'status': 'outside_cell_centres'}
    i = indices[0]
    x = centers[i]+slab.dx*(f[i]-.5)/(f[i]-f[i+1])
    return {'position_m': float(x), 'status': 'resolved'}


def trajectory(material, slab, initial, schedule, *, method, backend=None,
               sample_steps=None):
    """Integrate end-of-step wall/source inputs on one explicit common schedule.

    Each schedule entry is (dt, left, right, source). Stored energy is also
    recomputed from initial/final phase fractions, independently of summed
    step-ledger storage. Positive boundary/source heat enters the slab.
    """
    solver = METHODS[method]
    old = np.full(slab.cells, initial) if np.ndim(initial) == 0 else np.array(initial, copy=True)
    initial_state = old.copy()
    initial_fraction = material.liquid_fraction(initial_state)
    total_scale = slab.area*slab.length*(material.capacity*(material.liquidus-material.solidus)
                                       + material.latent_density)
    heat, defects, elapsed = [], [], []
    samples = []
    max_residual = max_defect = max_cumulative = max_closure = 0.
    for index, (dt, left, right, source) in enumerate(schedule, 1):
        step = solver(material, slab, old, dt, source=source, left=left, right=right, backend=backend)
        heat.append(fsum((step.ledger.left_J, step.ledger.right_J, step.ledger.source_J)))
        defects.append(step.ledger.defect_J)
        elapsed.append(dt)
        stored = slab.area*slab.dx*fsum(
            material.capacity*(np.asarray(step.temperature)-initial_state)
            + material.latent_density*(np.asarray(step.liquid_fraction)-initial_fraction))
        cumulative = fsum((stored, -fsum(heat)))
        closure = abs(cumulative-fsum(defects))/total_scale
        max_residual = max(max_residual, step.scaled_residual)
        max_defect = max(max_defect, step.ledger.scaled_defect)
        max_cumulative = max(max_cumulative, abs(cumulative)/total_scale)
        max_closure = max(max_closure, closure)
        if sample_steps is None or index in sample_steps:
            samples.append({'time_s': fsum(elapsed), **asdict(step),
                            'cumulative_heat_J': fsum(heat),
                            'cumulative_defect_J': cumulative,
                            'front': front_position(step.liquid_fraction, slab)})
        old = step.temperature
    if not samples:
        raise ValueError('schedule must produce at least one sample')
    summary = dict(max_scaled_equation_residual=max_residual,
                   max_scaled_step_energy_defect=max_defect,
                   max_scaled_cumulative_energy_defect=max_cumulative,
                   max_scaled_ledger_closure=max_closure)
    if not all(isfinite(v) for v in summary.values()) or max_residual > 1e-10 or max_closure > 1e-10:
        raise RuntimeError('comparison implementation/ledger gate failed')
    if method == 'enthalpy' and max_defect > 1e-9:
        raise RuntimeError('conservative energy gate failed')
    return {'samples': samples, 'steps': len(elapsed), **summary}


def experiment(backend=None):
    config = json.loads((ROOT/'phasechange/benchmarks.json').read_text())
    protocol = json.loads((ROOT/'phasechange/comparison.json').read_text())
    material = Material(**config['material'])
    params = {k: config['material'][k] for k in
              ('density', 'heat_capacity', 'latent_heat', 'solidus', 'liquidus')}
    uniform = []
    for case in config['uniform_cases']:
        for subdivisions in config['subdivisions']:
            schedule, reference = [], []
            start = case['initial']
            for duration, power in zip(case['durations'], case['powers']):
                schedule.extend([(duration/subdivisions, None, None, power)]*subdivisions)
                reference.extend(event_temperature(start, duration*j/subdivisions, power, **params)
                                 for j in range(1, subdivisions+1))
                start = reference[-1]
            row = {'case': case['name'], 'subdivisions': subdivisions,
                   'reference_temperature_K': reference, 'methods': {}}
            for name in METHODS:
                result = trajectory(material, Slab(1., 1), case['initial'], schedule,
                                    method=name, backend=backend)
                temperatures = [s['temperature'][0] for s in result['samples']]
                result['max_temperature_error_K'] = float(np.max(np.abs(np.asarray(temperatures)-reference)))
                result['final_temperature_error_K'] = abs(temperatures[-1]-reference[-1])
                if name == 'enthalpy' and result['max_temperature_error_K'] > 1e-9:
                    raise RuntimeError('uniform independent reference failed')
                row['methods'][name] = result
            uniform.append(row)

    settings = protocol['heat_cool_slab']
    slab = Slab(**{k: settings[k] for k in ('length', 'cells', 'area')})
    slab_cases = []
    for subdivisions in settings['subdivisions']:
        schedule = [(settings['base_dt_s']/subdivisions, segment['left_K'], None, segment['source_W_m3'])
                    for segment in settings['segments']
                    for _ in range(settings['steps_per_segment']*subdivisions)]
        methods = {name: trajectory(material, slab, settings['initial_K'], schedule, method=name,
                                   backend=backend, sample_steps=set(range(subdivisions, len(schedule)+1, subdivisions)))
                   for name in METHODS}
        differences = [float(np.max(np.abs(np.asarray(a['temperature'])-b['temperature'])))
                       for a,b in zip(methods['enthalpy']['samples'], methods['frozen']['samples'])]
        slab_cases.append({'subdivisions': subdivisions, 'methods': methods,
                           'max_sampled_method_difference_K': max(differences)})

    plan = config['stefan_plan']; front_settings = protocol['front_comparison']
    front_material = replace(material, solidus=plan['melting_temperature_K'],
                             liquidus=plan['melting_temperature_K']+front_settings['width_K'])
    front_slab = Slab(plan['domain_length_m'], front_settings['cells'], plan['area_m2'])
    front_cases = []
    for dt in front_settings['dt_s']:
        steps = [round(t/dt) for t in plan['observation_times_s']]
        if any(abs(i*dt-t) > 1e-12 for i,t in zip(steps, plan['observation_times_s'])):
            raise ValueError('observation times must lie on the shared time grid')
        schedule = [(dt, plan['hot_wall_K'], None, 0.)]*steps[-1]
        methods = {name: trajectory(front_material, front_slab, plan['melting_temperature_K'],
                                   schedule, method=name, backend=backend, sample_steps=set(steps))
                   for name in METHODS}
        front_cases.append({'dt_s': dt, 'methods': methods})

    files = sorted((ROOT/'phasechange').glob('*.py'))+sorted((ROOT/'phasechange').glob('*.json'))
    files += sorted(p for p in (ROOT/'native').glob('*') if p.is_file())
    return {'schema_version': 1,
            'scope': 'Matched-method verification and coarse front integration; front differences are not errors against an independent reference. No continuum, physical-validity or performance claim.',
            'backend': 'python' if backend is None else 'native',
            'frozen_endpoint_rule': 'right-sided tangent fixed at the old state, for heating and cooling',
            'energy_normalization': 'slab volume times [rho*c*(liquidus-solidus)+rho*L]; maximum cumulative defect covers all steps, including unsampled steps',
            'source_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'protocol': protocol, 'benchmark': config,
            'uniform': uniform, 'heat_cool_slab': slab_cases, 'front_comparison': front_cases}
