# Enthalpy conservation across phase endpoints

## Research contract

When a thermal step crosses a melting interval, how do conservative enthalpy
increments and frozen apparent heat capacity differ in energy balance,
temperature and front accuracy? This is a verification and synthesis study,
not a claim of a new discretization or a calibrated material model.

The November 3-4 requested break is excluded. Active day 1 starts November 5,
2026; target release November 11, earliest November 9, justified repairs only
November 12-13. A failed day-9 gate retains this project with a recovery plan.

## Model, units and equations

Constant density rho [kg/m3], specific heat c [J/kg/K], conductivity k [W/m/K],
latent heat L [J/kg], temperatures in kelvin. With solidus Ts < Tl liquidus,

```text
f(T) = 0                      T <= Ts
       (T-Ts)/(Tl-Ts)          Ts < T < Tl
       1                      T >= Tl
H(T) = rho*c*(T-Tref) + rho*L*f(T)             [J/m3]
dH/dt = div(k grad T) + q                     [W/m3]
C_m = rho*c + rho*L/(Tl-Ts)                   [J/m3/K]
```

H is continuous and strictly increasing. Its inverse has solid, mushy and
liquid affine branches. `Material.tangent(T, side=...)` chooses a one-sided
derivative at kinks; no classical two-sided derivative exists there. L=0 is
valid. A zero-width transition is excluded from this constitutive API.
Negative kelvin states, nonfinite inputs and overflowing outputs are errors.
Phase saturation is the constitutive definition, not clipping of an invalid
solver state. Tref changes stored enthalpy by a constant, never heat differences.

For a cell of volume V_i, the planned backward-Euler residual is

```text
R_i = V_i [H(T_i^{n+1})-H(T_i^n)]/dt
      - sum_j G_ij (T_j^{n+1}-T_i^{n+1}) - Q_boundary,i - V_i q_i.
```

G_ij is a positive face conductance [W/K]. The independent global ledger
recomputes sum_i V_i delta H_i minus dt times boundary/source power; interior
fluxes cancel. The planned comparator replaces delta H by
C(T_i^n)(T_i^{n+1}-T_i^n), with the same mesh, boundary inputs and time grid.
An unfavorable result or a tie is retained.

## Architecture and current implementation

Original `phasechange/` code is separate from six unmodified versioned module
snapshots. `Material` is immutable; vectorized constitutive operations use
NumPy. `integrate_uniform` returns an immutable `UniformHistory` with time,
temperature, enthalpy, phase fraction, supplied heat and recomputed energy
defect. The exact update H_new=H_old+q*dt handles multiple phase crossings.
The independent scalar oracle spends time to reach each phase boundary, then
changes dT/dt; it calls none of the constitutive methods or integrator.

Next: immutable finite-volume mesh, safeguarded nonlinear solver and shared
step result; then original C++17 residual/Jacobian/tridiagonal kernels with an
explicit Python interface. C++ and spatial phase-change methods are pending.
The existing native component checks are regression evidence only.

```sh
python -m phasechange --output results/enthalpy-baseline.json
python -m pytest tests/test_phasechange.py -q
```

Run these from the repository root after `python suite.py install`. The CLI
checks numerical gates and writes deterministic JSON bound to SHA-256 hashes
of every phase-change source and the benchmark file. CI reproduces the JSON
on Linux and Windows; the regression test rejects stale source hashes.

## Fixed synthetic experiments and scales

`phasechange/benchmarks.json` is the versioned parameter authority. The uniform
baseline uses rho=1200, c=900, k=1.8, L=90000, Ts=309, Tl=311, Tref=300 in SI.
Six schedules cover single-phase heating, mushy heating/cooling, both endpoint
crossings in each direction and a closed heating/cooling cycle. Each is sampled
with 1, 2, 4, 8 and 16 subdivisions. Error scales are Tl-Ts=2 K and
rho*L+rho*c*(Tl-Ts)=110160000 J/m3. These are implementation and constitutive
checks: exact constant-power integration has no time-discretization error.

The future Stefan test is one-phase melting from 310 K with a 330 K left wall,
0.1 m length, unit area and an insulated far end. For the sharp semi-infinite
reference, alpha=k/(rho*c), Ste=c*(330-310)/L=0.2,

```text
sqrt(pi)*lambda*exp(lambda^2)*erf(lambda) = Ste
s(t) = 2*lambda*sqrt(alpha*t)
T(x,t) = 330 - 20*erf(x/(2*sqrt(alpha*t)))/erf(lambda),  x <= s(t)
T(x,t) = 310,                                         x > s(t).
```

The regularized test uses Ts=310 and Tl=310+width, so initial solid has zero
liquid fraction. The positive-time observation window is fixed at 30, 60, 90,
120 seconds. Independent mesh, time and width sweeps and doubled-domain checks
are specified in JSON before implementation. The f=0.5 front is interpolated
between cell centers. This defines a regularized front whose error against a
sharp interface includes model regularization as well as discrete error.

## Acceptance gates for the complete cycle

1. Constitutive inverse and independent uniform heat/cool errors <=1e-12 on
   declared scales; monotonicity, endpoints and fractions in [0,1].
2. Converged step scaled residual <=1e-10 and independently recomputed global
   energy defect <=1e-9; iteration failure explicit.
3. Python/C++ residual/Jacobian agreement <=1e-11, state agreement <=1e-9;
   finite-difference checks away from kinks, directional endpoint and native
   invariant checks. Scaling must be documented before applying a threshold.
4. Zero-latent limit agrees with independent heat references; four-level time
   order 0.9-1.1 and smooth spatial order 1.8-2.2. Exact semidiscrete evolution
   separates temporal error from continuum spatial error.
5. Stefan front error <=2% at all four observation times at final settings;
   separate mesh/time/width/domain effects. Investigate failure rather than
   selecting favorable times or silently retuning parameters.
6. Matched frozen-capacity comparison records energy, temperature and front
   errors, including endpoint crossings; no universal-superiority claim.
7. Existing regressions and all snapshot hashes pass; Linux/Windows reproduce
   connected Python/C++ experiments with source-bound evidence.
8. An original synthesis at `paper/paper.tex` and `paper/paper.pdf`, <=5 pages
   including references, reports measured results. Compile, check diagnostics,
   render and inspect every page; manuscript gates in CI. The existing v0.1.0
   paper remains supporting prior work until the replacement is ready.
9. Merge/release only after all gates pass. No physical validation, performance
   or material-calibration claims without corresponding independent evidence.

## Provenance, limitations and references

DTU MMME (Modelling in Materials and Manufacturing Engineering) solidification,
phase-fraction and implicit-conduction topics motivate the question. Newly
available solidification/latent-heat and implicit-method pages were inspected
read-only on November 5. Only topic provenance is used; no private text,
geometry, exercise parameters, code, solutions or data are included. The new
ATSA grey-box lecture is a candidate for the later Vibration cycle.

This model excludes fluid flow, segregation, nucleation, hysteresis, variable
properties, deformation and physical experimental validation. A finite melting
interval is an intentional regularization, not proof of alloy physics. A small
nonlinear residual alone establishes neither continuum nor physical accuracy.

Enthalpy formulations are established methods: V. R. Voller, M. Cross and
N. C. Markatos, "An enthalpy method for convection/diffusion phase change,"
International Journal for Numerical Methods in Engineering 24 (1987), 271-284,
https://doi.org/10.1002/nme.1620240119. Bibliographic details checked against the
publisher and author institution; no implementation is copied from this work.
