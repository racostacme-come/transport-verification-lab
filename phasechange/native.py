"""Explicit C++17 assembly/LDL backend via C ABI 1, without fallback."""
import ctypes
import os
from pathlib import Path
import numpy as np
from .finite_volume import _field


def _owned(values):
    return np.array(values, dtype=np.float64, order='C', copy=True)


def _check(status):
    if status == 1:
        raise ValueError('native kernel rejected invalid input')
    if status == 2:
        raise FloatingPointError('native arithmetic exceeds representable range')
    if status == 3:
        raise ValueError('native LDL requires a positive definite matrix')
    if status:
        raise RuntimeError(f'native kernel failed with status {status}')


class NativeKernel:
    """Own a library handle; each call copies input buffers and owns outputs.

    Load library=path, PHASECHANGE_LIBRARY, or suite.py's ignored build record.
    Native calls retain no pointers and use no mutable global state.
    """
    def __init__(self, library=None):
        path = library or os.environ.get('PHASECHANGE_LIBRARY')
        record = Path(__file__).resolve().parents[1]/'build/phasechange-library.txt'
        if not path and record.exists():
            path = record.read_text(encoding='utf-8').strip()
        if not path:
            raise RuntimeError('Run python suite.py native or set PHASECHANGE_LIBRARY')
        self.path = Path(path).resolve(strict=True)
        self._library = ctypes.CDLL(str(self.path))
        version = self._library.pc_abi_version
        version.argtypes = []; version.restype = ctypes.c_int
        if version() != 1:
            raise RuntimeError('Expected phase-change ABI 1')
        array = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags=('C_CONTIGUOUS','ALIGNED'))
        self._assemble = self._library.pc_assemble
        self._assemble.argtypes = ([ctypes.c_size_t]+[ctypes.c_double]*8+
            [ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_double,ctypes.c_int]+[array]*6)
        self._assemble.restype = ctypes.c_int
        self._solve = self._library.pc_solve
        self._solve.argtypes = [ctypes.c_size_t]+[array]*4
        self._solve.restype = ctypes.c_int

    def assemble(self, material, slab, previous, trial, dt, *, source=0., left=None, right=None, side='right'):
        if side not in ('left','right'):
            raise ValueError('side must be left or right')
        n = slab.cells
        old, temp, power = [_owned(_field(v,n,name)) for v,name in
                            ((previous,'previous'),(trial,'trial'),(source,'source'))]
        residual, diagonal, off = np.empty(n), np.empty(n), np.empty(n-1)
        _check(self._assemble(n,material.capacity,material.latent_density,material.solidus,
            material.liquidus,material.conductivity,slab.dx,slab.area,dt,
            int(left is not None),0. if left is None else left,
            int(right is not None),0. if right is None else right,int(side=='left'),
            old,temp,power,residual,diagonal,off))
        return residual,diagonal,off

    def solve(self, diagonal, off, rhs):
        d,e,b = map(_owned,(diagonal,off,rhs))
        if d.ndim!=1 or len(d)==0 or e.shape!=(len(d)-1,) or b.shape!=d.shape:
            raise ValueError('invalid tridiagonal shapes')
        result = np.empty_like(d)
        _check(self._solve(len(d),d,e,b,result))
        return result
