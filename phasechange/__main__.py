"""Run a source-bound baseline or matched-method comparison."""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from .material import Material
from .reference import event_temperature
from .uniform import integrate_uniform

ROOT = Path(__file__).resolve().parents[1]


def experiment():
    config = json.loads((ROOT/'phasechange/benchmarks.json').read_text(encoding='utf-8'))
    material = Material(**config['material'])
    params = {k: config['material'][k] for k in ('density', 'heat_capacity', 'latent_heat', 'solidus', 'liquidus')}
    tscale = material.liquidus-material.solidus
    hscale = material.latent_density+material.capacity*tscale
    rows = []
    for case in config['uniform_cases']:
        for subdivisions in config['subdivisions']:
            dt = [d/subdivisions for d in case['durations'] for _ in range(subdivisions)]
            powers = [q for q in case['powers'] for _ in range(subdivisions)]
            history = integrate_uniform(material, case['initial'], dt, powers)
            # Oracle evaluates from each original segment start, not recursively
            # through numerical samples; no shared inverse or cumulative error.
            reference = [case['initial']]
            start = case['initial']
            for duration, power in zip(case['durations'], case['powers']):
                reference.extend(event_temperature(start, duration*j/subdivisions, power, **params)
                                 for j in range(1, subdivisions+1))
                start = reference[-1]
            error = float(np.max(np.abs(np.asarray(history.temperature)-reference))/tscale)
            energy = float(np.max(np.abs(history.energy_defect))/hscale)
            rows.append({'case': case['name'], 'subdivisions': subdivisions,
                         'history': asdict(history), 'reference_temperature_K': reference,
                         'scaled_temperature_error': error, 'scaled_energy_defect': energy})
    files = [ROOT/'phasechange'/name for name in ('__init__.py', '__main__.py',
             'material.py', 'reference.py', 'uniform.py', 'benchmarks.json')]
    return {'schema_version': 1, 'scope': 'Uniform constitutive verification only; no spatial solver, front accuracy or physical validation claim.',
            'source_sha256': {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'temperature_scale_K': tscale, 'enthalpy_scale_J_m3': hscale,
            'cases': rows, 'max_scaled_temperature_error': max(r['scaled_temperature_error'] for r in rows),
            'max_scaled_energy_defect': max(r['scaled_energy_defect'] for r in rows)}


def check(data):
    for key in ('max_scaled_temperature_error', 'max_scaled_energy_defect'):
        if not np.isfinite(data[key]) or data[key] > 1e-12:
            raise RuntimeError(f"uniform acceptance failed: {key}={data[key]}")
    for row in data['cases']:
        if not all(0 <= f <= 1 for f in row['history']['liquid_fraction']):
            raise RuntimeError('phase fraction out of range')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', choices=('baseline', 'comparison'), default='baseline')
    parser.add_argument('--backend', choices=('python', 'native'), default='python')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--check-against', type=Path, help='Comparison JSON to verify before writing')
    args = parser.parse_args()
    if args.study == 'baseline':
        if args.check_against:
            parser.error('--check-against applies to --study comparison')
        if args.backend != 'python':
            parser.error('baseline uses the exact Python constitutive inverse; choose --study comparison for native')
        data = experiment(); check(data)
        summary = {k: v for k, v in data.items() if k not in ('source_sha256', 'cases')}
    else:
        from .comparison import experiment as compare
        from .native import NativeKernel
        data = compare(None if args.backend == 'python' else NativeKernel())
        if args.check_against:
            from .comparison import verify_reproduction
            verify_reproduction(data, json.loads(args.check_against.read_text(encoding='utf-8')))
        summary = {'study': args.study, 'backend': args.backend,
                   'uniform_cases': len(data['uniform']), 'slab_cases': len(data['heat_cool_slab']),
                   'front_cases': len(data['front_comparison'])}
    suffix = '-native' if args.backend == 'native' else ''
    output = args.output or ROOT/f'results/enthalpy-{args.study}{suffix}.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
