"""Reproduce native assembly and connected trajectory checks (no speed claims)."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import numpy as np
from .finite_volume import Slab, advance
from .material import Material
from .native import NativeKernel

ROOT = Path(__file__).resolve().parents[1]


def reference_assembly(mat, slab, previous, trial, dt, source, left, right, side='right'):
    """Scalar cell-balance oracle; integrate capacities between old/trial T.

    Does not call native assembly or Material.enthalpy/fraction/tangent.
    Neighbor contributions are assembled per cell, not shared face updates.
    """
    n = slab.cells; volume = slab.area*slab.dx
    conductance = mat.conductivity*slab.area/slab.dx
    residual = np.empty(n); diagonal = np.empty(n)
    for i, (old, t, q) in enumerate(zip(previous, trial, source)):
        low, high = sorted((old,t))
        mushy_span = max(0., min(high,mat.liquidus)-max(low,mat.solidus))
        sign = 1. if t>=old else -1.
        increment = sign*(mat.capacity*(high-low)+mat.latent_density*mushy_span/(mat.liquidus-mat.solidus))
        heat_in = 0.; neighbor_g = 0.
        for j in (i-1,i+1):
            if 0<=j<n:
                heat_in += conductance*(trial[j]-t); neighbor_g += conductance
        for active,wall in ((i==0,left),(i==n-1,right)):
            if active and wall is not None:
                heat_in += 2*conductance*(wall-t); neighbor_g += 2*conductance
        inside = (mat.solidus<=t<mat.liquidus) if side=='right' else (mat.solidus<t<=mat.liquidus)
        tangent = mat.capacity+(mat.latent_density/(mat.liquidus-mat.solidus) if inside else 0.)
        residual[i] = volume*increment-dt*(heat_in+volume*q)
        diagonal[i] = volume*tangent+dt*neighbor_g
    return residual,diagonal,np.full(n-1,-dt*conductance)


def experiment(kernel=None):
    kernel = kernel or NativeKernel()
    material = Material(**json.loads((ROOT/'phasechange/benchmarks.json').read_text())['material'])
    rng = np.random.default_rng(712023)
    rows = []
    for n in (1,2,24):
        for latent in (0., material.latent_heat):
            mat = replace(material, latent_heat=latent); slab = Slab(.012,n,.003)
            for left,right in ((None,None),(335.,None),(335.,295.)):
                for dt in (.25,10.):
                    old = rng.uniform(300,320,n); trial = rng.uniform(300,320,n)
                    source = rng.uniform(-2.7e6,2.7e6,n)
                    actual = kernel.assemble(mat,slab,old,trial,dt,source=source,left=left,right=right)
                    reference = reference_assembly(mat,slab,old,trial,dt,source,left,right)
                    energy_scale = slab.area*slab.dx*(mat.capacity*(mat.liquidus-mat.solidus)+mat.latent_density)
                    jacobian_scale = slab.area*slab.dx*mat.mushy_capacity+4*dt*mat.conductivity*slab.area/slab.dx
                    rerror = float(np.max(np.abs(actual[0]-reference[0]))/energy_scale)
                    jerror = float(max(np.max(np.abs(actual[1]-reference[1])),
                                      np.max(np.abs(actual[2]-reference[2]),initial=0))/jacobian_scale)
                    py = advance(mat,slab,old,dt,source=source,left=left,right=right)
                    native = advance(mat,slab,old,dt,source=source,left=left,right=right,backend=kernel)
                    terror = float(np.max(np.abs(np.asarray(py.temperature)-native.temperature))/(mat.liquidus-mat.solidus))
                    rows.append(dict(cells=n,latent_heat=latent,left_K=left,right_K=right,dt_s=dt,
                                     scaled_residual_difference=rerror,scaled_jacobian_difference=jerror,
                                     scaled_state_difference=terror,scaled_energy_defect=native.ledger.scaled_defect))
    slab = Slab(.012,24,.003); py_state=np.full(24,305.); cpp_state=py_state.copy(); trajectory=[]
    for wall,power in ((335.,2.7e6),(295.,-2.7e6)):
        for _ in range(5):
            py = advance(material,slab,py_state,10.,left=wall,source=power)
            cpp = advance(material,slab,cpp_state,10.,left=wall,source=power,backend=kernel)
            trajectory.append({'python':asdict(py),'native':asdict(cpp),
                'scaled_state_difference':float(np.max(np.abs(np.asarray(py.temperature)-cpp.temperature))/2.)})
            py_state=py.temperature; cpp_state=cpp.temperature
    paths = [ROOT/'phasechange'/x for x in ('native.py','native_checks.py','finite_volume.py','material.py','benchmarks.json')]
    paths += sorted((ROOT/'native').glob('*'))
    data={'schema_version':1,'scope':'Native implementation agreement, independent assembly oracle and shared Python safeguard/ledger; no independent nonlinear algorithm, continuum or performance claim.',
          'source_sha256':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()},
          'assembly_and_state_cases':rows,'connected_trajectory':trajectory,
          'max_scaled_residual_difference':max(x['scaled_residual_difference'] for x in rows),
          'max_scaled_jacobian_difference':max(x['scaled_jacobian_difference'] for x in rows),
          'max_scaled_state_difference':max(x['scaled_state_difference'] for x in rows+trajectory),
          'max_scaled_energy_defect':max([x['scaled_energy_defect'] for x in rows]+[x['native']['ledger']['scaled_defect'] for x in trajectory])}
    for key,limit in (('max_scaled_residual_difference',1e-11),('max_scaled_jacobian_difference',1e-11),
                      ('max_scaled_state_difference',1e-9),('max_scaled_energy_defect',1e-9)):
        if not np.isfinite(data[key]) or data[key]>limit:
            raise RuntimeError(f'native gate failed: {key}={data[key]}')
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/enthalpy-native.json')
    args=parser.parse_args(); data=experiment()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print({k:v for k,v in data.items() if k.startswith('max_')})
