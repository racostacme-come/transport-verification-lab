# Matched enthalpy and frozen-capacity comparison

The conservative solver uses the full enthalpy increment. The deliberately
approximate comparator freezes apparent capacity at the previous temperature:

```text
V C(T_old) (T_new - T_old) = dt [sum_faces G (T_neighbor,new-T_new) + V q]
C(T) = rho c + rho L/(Tl-Ts)  for Ts <= T < Tl
       rho c                 otherwise
```

The endpoint convention is the right-sided tangent for both heating and
cooling. It is explicit and fixed, not a direction-dependent correction.
For example, cooling from exactly Ts uses the mushy capacity for that whole
step; cooling from exactly Tl uses the sensible capacity. This asymmetry can
affect errors. No smoothing, endpoint substeps or energy correction are hidden
inside this comparator. It represents this particular frozen algorithm,
not all apparent-capacity methods.

The same `Material`, `Slab`, source, wall temperatures, time steps and `Step`
result contract serve both methods. `advance_frozen` solves one positive
definite tridiagonal system in temperature increments. In native mode the
existing original C++17 kernel supplies the Jacobian at the old state and the
LDL solve. Python reconstructs the accepted linear balance independently.
Both backends reject invalid/nonfinite inputs and negative Kelvin results.
There is no fallback from a missing native library.

```python
from phasechange import Material
from phasechange.finite_volume import Slab, advance
from phasechange.frozen import advance_frozen
from phasechange.native import NativeKernel

mat = Material(1200, 900, 1.8, 90000, 309, 311, reference=300)
mesh = Slab(.012, 24, .003)
conservative = advance(mat, mesh, 305., 10., left=335., source=2.7e6)
frozen = advance_frozen(mat, mesh, 305., 10., left=335., source=2.7e6,
                        backend=NativeKernel())
```

`Step.scaled_residual` is the equation residual for the selected method. For
frozen capacity it is **not** the enthalpy residual. Both `Step.ledger` objects
independently integrate the actual sensible/latent heat between accepted
states and reconstruct wall/source heat. A large frozen energy defect is a
measured outcome, not a solver failure. Both algebraic residuals must be at
most 1e-10 on the cell energy scale; only the conservative method has the
1e-9 true-energy gate.

`comparison.trajectory` accumulates signed wall/source inputs and recomputes
the total stored-energy change from initial and current temperature/fraction.
It checks this against the sum of independent step defects. The normalization
is total volume times `rho*c*(Tl-Ts)+rho*L`. It is an energy scale, not the
relative error against supplied heat, which can vanish in a cycle. Maximum
cumulative defect includes every step, even when profiles are sparsely sampled.

## Reproduction

After `python suite.py install` and `python suite.py native`:

```sh
python -m phasechange --study comparison --output results/enthalpy-comparison.json
python -m phasechange --study comparison --backend native --check-against results/enthalpy-comparison.json --output validation/comparison-native.json
python -m pytest tests/test_frozen_comparison.py -q
```

The original `python -m phasechange` uniform baseline remains available.
The new CLI rejects native mode for that exact constitutive baseline, rather
than implying it uses a C++ time integrator. Native comparison defaults to a
separate output filename to preserve the Python reference.

The day-1 `benchmarks.json` is unchanged. `comparison.json` fixes the day-4
protocol: 30 uniform schedule/subdivision combinations; five 24-cell heating/
cooling runs at matched sample times; four coarse front runs at N=100, width
1 K and dt=0.5, 0.25, 0.125, 0.0625 s. The latter use the day-1 wall, initial
temperature, domain and four observation times. Every row retains both
methods. Uniform errors use the independent phase-event reference. Slab and
front results are method differences, **not independent accuracy estimates**.

The front is the unique decreasing liquid-fraction 0.5 crossing interpolated
between cell centres. A 0.5 plateau, nonmonotone field or crossing outside the
centre interval gets an explicit status and null position, without wall
extrapolation. The JSON includes sampled profiles and all-step energy maxima,
parameter copies and SHA-256 hashes of all original Python/native study files.
CI reproduces Python histories against the versioned evidence, then native
histories against Python. It checks temperature within 2e-9 K, fraction/front
position within 1e-9, and dimensional energies within 2e-7 absolute plus 2e-9
relative. Residual/closure gates are checked independently of platform
roundoff and iteration counts. These are implementation reproduction gates.

Independent sharp Stefan comparison and separated mesh/time/width/domain
studies remain day-5 work. The coarse fixed-grid sweep does not replace those
acceptance gates. No physical validation or universal method superiority is
claimed. Topic provenance and model exclusions remain in the
[research contract](enthalpy-conservation.md).
