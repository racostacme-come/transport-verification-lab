"""Independent heat and fixed day-1 Stefan studies, on either kernel backend.

Run python -m phasechange.refinement. Failed accuracy gates are written to the
output before a nonzero exit; failed nonlinear solves are never accepted.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np

from .comparison import trajectory, verify_reproduction
from .continuum import cosine_reference, similarity_parameter, stefan_solution
from .finite_volume import Slab, advance
from .material import Material

ROOT = Path(__file__).resolve().parents[1]


def heat_study(material, backend=None):
    mat = replace(material, latent_heat=0.)
    alpha = mat.conductivity/mat.capacity
    length = .1
    final = length**2/(alpha*np.pi**2)  # one continuum decay time
    temporal = []
    mesh = Slab(length, 32)
    initial = 310+8*cosine_reference(mesh.cells, length, alpha, 0.)
    exact = 310+8*cosine_reference(mesh.cells, length, alpha, final, discrete=True)
    rate = 4*alpha*(np.sin(np.pi/(2*mesh.cells))/mesh.dx)**2
    for steps in (20, 40, 80, 160):
        state = initial.copy()
        for _ in range(steps):
            state = advance(mat, mesh, state, final/steps, backend=backend).temperature
        algebraic = 310+(initial-310)/(1+final/steps*rate)**steps
        temporal.append({'steps': steps, 'dt_s': final/steps,
                         'error_K': float(np.max(np.abs(np.asarray(state)-exact))),
                         'algebraic_error_K': float(np.max(np.abs(np.asarray(state)-algebraic)))})
    spatial = []
    for cells in (10, 20, 40, 80):
        semi = cosine_reference(cells, length, alpha, final, discrete=True)
        exact = cosine_reference(cells, length, alpha, final)
        spatial.append({'cells': cells, 'dx_m': length/cells,
                        'error_K': float(8*np.max(np.abs(semi-exact)))})
    for rows in (temporal, spatial):
        for coarse, fine in zip(rows, rows[1:]):
            fine['order'] = float(np.log2(coarse['error_K']/fine['error_K']))
    return {'length_m': length, 'final_time_s': final, 'amplitude_K': 8.,
            'temporal_cells': 32, 'temporal': temporal, 'spatial': spatial,
            'spatial_policy': 'Exact semidiscrete evolution versus continuum cell averages; no time error.'}


def front_case(material, plan, cells, dt, width, length, method, backend):
    mat = replace(material, solidus=plan['melting_temperature_K'],
                  liquidus=plan['melting_temperature_K']+width)
    slab = Slab(length, cells, plan['area_m2'])
    steps = [round(t/dt) for t in plan['observation_times_s']]
    if any(abs(i*dt-t) > 1e-12 for i,t in zip(steps, plan['observation_times_s'])):
        raise ValueError('observation times must lie on the fixed time grid')
    result = trajectory(mat, slab, mat.solidus,
                        [(dt,plan['hot_wall_K'],None,0.)]*steps[-1],
                        method=method, backend=backend, sample_steps=set(steps))
    ste = mat.heat_capacity*(plan['hot_wall_K']-mat.solidus)/mat.latent_heat
    samples = []
    for sample in result.pop('samples'):
        sharp, _ = stefan_solution([], sample['time_s'],
            diffusivity=mat.conductivity/mat.capacity, stefan_number=ste,
            melting=mat.solidus, wall=plan['hot_wall_K'])
        position = sample['front']['position_m']
        samples.append({'time_s':sample['time_s'], 'front':sample['front'],
            'sharp_front_m':sharp,
            'relative_front_error':None if position is None else abs(position/sharp-1),
            'temperature':sample['temperature'], 'liquid_fraction':sample['liquid_fraction'],
            'cumulative_defect_J':sample['cumulative_defect_J']})
    return {'cells':cells,'dt_s':dt,'width_K':width,'length_m':length,
            'method':method,'samples':samples,**result}


def acceptance(data):
    heat = data['heat']
    refined = data['cases']['refined']['enthalpy']
    errors = [s['relative_front_error'] for s in refined['samples']]
    return {
        'heat_temporal_order':all(.9 <= r['order'] <= 1.1 for r in heat['temporal'][1:]),
        'heat_spatial_order':all(1.8 <= r['order'] <= 2.2 for r in heat['spatial'][1:]),
        'heat_algebraic':all(r['algebraic_error_K'] <= 1e-9 for r in heat['temporal']),
        'refined_stefan_front':all(e is not None and e <= .02 for e in errors),
        'conservative_cumulative_energy':all(
            row['enthalpy']['max_scaled_cumulative_energy_defect'] <= 1e-9
            for row in data['cases'].values()),
    }


def experiment(backend=None):
    config = json.loads((ROOT/'phasechange/benchmarks.json').read_text())
    plan = config['stefan_plan']
    material = Material(**config['material'])
    n,dt,width,length = (plan['mesh_cells'][-1], plan['time_steps_s'][-1],
                         plan['melting_widths_K'][-1], plan['domain_length_m'])
    specifications = {'refined':(n,dt,width,length)}
    groups = {}
    for group, values in (('mesh',plan['mesh_cells']),('time',plan['time_steps_s']),
                          ('width',plan['melting_widths_K'])):
        names = []
        for i,value in enumerate(values):
            name = 'refined' if i == len(values)-1 else f'{group}-{i}'
            names.append(name)
            if name != 'refined':
                specifications[name] = (value if group=='mesh' else n,
                    value if group=='time' else dt, value if group=='width' else width, length)
        groups[group] = names
    specifications['double-domain'] = (2*n,dt,width,2*length)
    cases = {}
    for name, settings in specifications.items():
        cases[name] = {method:front_case(material,plan,*settings,method,backend)
                       for method in ('enthalpy','frozen')}
        print(f"completed {name}", flush=True)
    domain = {}
    for method in ('enthalpy','frozen'):
        base, extended = cases['refined'][method],cases['double-domain'][method]
        pairs = list(zip(base['samples'],extended['samples']))
        resolved = all(a['front']['position_m'] is not None and b['front']['position_m'] is not None
                       for a,b in pairs)
        domain[method] = {
            'front_status':'resolved' if resolved else 'unresolved',
            'max_front_difference_m':max(abs(a['front']['position_m']-b['front']['position_m'])
                for a,b in pairs) if resolved else None,
            'max_common_temperature_difference_K':max(float(np.max(np.abs(
                np.asarray(a['temperature'])-np.asarray(b['temperature'])[:n])))
                for a,b in zip(base['samples'],extended['samples']))}
    paths = sorted((ROOT/'phasechange').glob('*.py'))+sorted((ROOT/'phasechange').glob('*.json'))
    paths += sorted(p for p in (ROOT/'native').glob('*') if p.is_file())
    data = {'schema_version':1,'backend':'python' if backend is None else 'native',
        'scope':'Numerical verification against independent analytical references; no physical validation.',
        'source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'protocol':plan,'groups':groups,'cases':cases,'domain_sensitivity':domain,
        'similarity_parameter':similarity_parameter(material.heat_capacity*(plan['hot_wall_K']-
            plan['melting_temperature_K'])/material.latent_heat),
        'heat':heat_study(material,backend)}
    data['acceptance'] = acceptance(data)
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend',choices=('python','native'),default='python')
    parser.add_argument('--output',type=Path,default=ROOT/'results/enthalpy-refinement.json')
    parser.add_argument('--check-against',type=Path)
    args = parser.parse_args()
    from .native import NativeKernel
    data = experiment(None if args.backend=='python' else NativeKernel())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    if args.check_against:
        verify_reproduction(data,json.loads(args.check_against.read_text()))
    print(json.dumps(data['acceptance']))
    if not all(data['acceptance'].values()):
        raise SystemExit('Accuracy gate failed; results retained for investigation.')


if __name__ == '__main__':
    main()
