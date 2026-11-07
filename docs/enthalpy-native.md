# Native enthalpy backend

Original C++17 code in `native/` is MIT licensed with the repository. It changes
no imported snapshots. Python owns the nonlinear iteration, line search,
material contract, acceptance tests and independent signed energy ledger.
The optional backend assembles the residual/Jacobian and solves its symmetric
positive definite tridiagonal system. Both paths use the same `Step` result.

```sh
python suite.py install
python suite.py native --build-root ../../.native-transport
python suite.py test
python -m phasechange.native_checks --output results/enthalpy-native.json
```

The native build records its absolute library path in ignored
`build/phasechange-library.txt`. CI builds before tests on Linux and Windows.
It runs a standalone native invariant executable plus Python binding tests;
a missing library fails the native tests instead of silently skipping them.
For a direct build, use `cmake -S native -B build/phasechange`,
`cmake --build build/phasechange --config Release`, then
`ctest --test-dir build/phasechange -C Release --output-on-failure`.
The generated `library-Release.txt` gives the library path for explicit loading.

```python
from phasechange import Material
from phasechange.finite_volume import Slab, advance
from phasechange.native import NativeKernel

material = Material(1200, 900, 1.8, 90000, 309, 311, reference=300)
slab = Slab(.012, 24, .003)
kernel = NativeKernel()  # or library=path, or PHASECHANGE_LIBRARY
step = advance(material, slab, 305., 10., left=335., source=2.7e6,
               backend=kernel)
r, diagonal, off = kernel.assemble(material, slab, 305., step.temperature,
                                   10., left=335., source=2.7e6)
correction = kernel.solve(diagonal, off, -r)
```

## Units, Jacobian and endpoint convention

The implemented residual is dt times the power residual in the research
contract: each entry is joules, and Jacobian entries are joules/kelvin.
For uniform volume V=A dx and g=dt k A/dx,

```text
R_i = V [rho*c*(T_i-old_i) + rho*L*(f(T_i)-f(old_i)) - dt*q_i]
      + sum_neighbor g*(T_i-T_neighbor) + wall terms.
J_ii = V*C(T_i) + sum_neighbor g + sum_present_wall 2g
J_i,i+1 = J_i+1,i = -g.
```

An absent wall is insulated. A prescribed wall adds 2g*(T_i-T_wall).
Both walls act on the same cell when n=1. At a phase kink `side='right'`
(default) means Ts is mushy and Tl is liquid; `side='left'` means Ts is
solid and Tl is mushy. The native Newton path uses the right derivative,
matching the Python path. Directional tests cover both conventions.

Negative finite trial temperatures extend the solid affine branch during line
search. Accepted states still pass the original nonnegative-kelvin validation;
no invalid final state is clipped. Old states must always be nonnegative.

The LDL^T solve requires positive finite pivots and has no pivoting or fallback.
For the intended Jacobian, positive storage plus positive conductances gives
positive definiteness. The general `solve` API rejects indefinite systems.
This is not a general-purpose nonsymmetric or indefinite tridiagonal solver.

## ABI and failure behavior

`native/kernels.h` defines C ABI version 1. Input/output buffers are disjoint,
caller-owned double arrays, length n except off-diagonal length n-1. No pointer
is retained. Python copies all inputs into owned contiguous float64 buffers.
C++ uses temporary storage and writes outputs only after a successful call.
No exception crosses the ABI. Status codes distinguish invalid input (1),
unrepresentable arithmetic (2), nonpositive LDL pivots (3), and internal or
allocation failure (4). ABI mismatch and missing libraries raise explicit errors.

## Reproducible verification and limits

The seeded 36-case study varies 1/2/24 cells, zero/nonzero latent heat,
insulated/mixed/two-wall boundaries and 0.25/10 second steps. An independent
scalar oracle integrates capacities over the traversed temperature interval
and assembles each cell balance separately; it shares no native constitutive
or face-update code. Converged states are compared with the Python default.
The ten-step day-2 trajectory is run independently through each backend.

Residual differences are normalized by V*(rho*c*(Tl-Ts)+rho*L); Jacobian
differences by V*C_m+4g; temperature differences by Tl-Ts. Gates are 1e-11,
1e-11 and 1e-9, respectively. The independent global energy ledger retains its
1e-9 gate. Source hashes include kernels, native tests, CMake, binding, solver,
material, benchmark configuration and experiment code. The earlier slab
manifest is explicitly refreshed for the backend hook with identical numerical
JSON values; no imported provenance/integrity manifest is changed.

These are implementation checks. Python/C++ agreement shares nonlinear
safeguarding, so it cannot prove nonlinear or physical correctness by itself.
Continuum/front studies and frozen-capacity comparison remain separate gates.
No speedup, sanitizer coverage, independent native nonlinear solver or physical
validation is claimed. The original three-page thermal paper remains unchanged
until the phase-change synthesis is supported by later experiments.
