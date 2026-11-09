# Independent heat and Stefan verification

`python -m phasechange.refinement` runs the fixed day-1 front protocol and
isolated zero-latent heat studies. It writes `results/enthalpy-refinement.json`,
including source SHA-256 hashes, all four observed states of each case, error
and energy histories, and explicit acceptance booleans. An accuracy failure
writes the evidence and exits unsuccessfully. Nonlinear failures raise errors.
Unresolved front crossings remain explicit, including in domain comparisons.

```sh
python suite.py native
python -m phasechange.refinement --check-against results/enthalpy-refinement.json --output validation/refinement-python.json
python -m phasechange.refinement --backend native --check-against validation/refinement-python.json --output validation/refinement-native.json
```

The two backends use the same Python nonlinear safeguard, front measurement
and energy ledger. Their agreement checks implementation, not independence
of the nonlinear algorithm. Both OS jobs reproduce both backends in CI.

## Analytical references

`phasechange/continuum.py` imports neither the material implementation nor
the numerical solvers. For constant properties, alpha=k/(rho*c), a hot wall
Tw, and solid initially at Tm, the one-phase sharp Stefan solution is

\[
s(t)=2\lambda\sqrt{\alpha t},\quad
\lambda e^{\lambda^2}\operatorname{erf}(\lambda)
=\frac{c(T_w-T_m)}{L\sqrt\pi},
\]
\[
T(x,t)=T_w-(T_w-T_m)\frac{\operatorname{erf}(x/(2\sqrt{\alpha t}))}
{\operatorname{erf}(\lambda)}\quad(0\le x<s),\qquad T=T_m\quad(x\ge s).
\]

The scalar constant is solved by bisection of a monotone logarithmic equation.
Independent tests check the heat PDE by centered differences, wall and interface
temperatures, the one-sided Stefan flux jump, square-root front scaling and
prescribed similarity constants. Extreme Stefan-number tests guard overflow
and underflow in the reference, rather than changing the fixed experiment.

For an insulated slab, cosine cell averages have continuum decay rate
alpha*(pi/length)^2. The cell-centred finite-volume rate is
4*alpha*sin(pi/(2*N))^2/dx^2. Exact exponential evolution with the latter
isolates spatial error; backward Euler compared against that exponential
isolates time error. An independently assembled incidence matrix and dense
matrix exponential check the semidiscrete oracle. Exact integration checks
the initial cell averages. No time error contaminates the spatial table.

## Fixed protocol and interpretation

The synthetic material and Stefan protocol are unchanged in
`phasechange/benchmarks.json`: rho=1200 kg/m3, c=900 J/(kg K), k=1.8 W/(m K),
L=90000 J/kg, Tm=310 K, Tw=330 K, area=1 m2, insulated far wall at 0.1 m.
The regularized material has solidus=Tm, liquidus=Tm+width and starts with
zero liquid fraction. Observations are 30, 60, 90 and 120 s; the measured
front is the interpolated decreasing f=0.5 crossing between cell centers.

The four-level mesh, time and width studies change only one parameter:
N=100/200/400/800, dt=0.5/0.25/0.125/0.0625 s, width=1/0.5/0.25/0.125 K.
Unvaried parameters use N=800, dt=0.0625 s, width=0.125 K. A final domain
check uses 0.2 m and N=1600, preserving dx. Shared finest cases are computed
once: 11 distinct settings, each run with enthalpy and frozen capacity.

At fixed nonzero width, errors against the sharp solution include
regularization, discretization and front interpolation effects. These are
controlled sensitivity studies, not a mathematically complete additive error
decomposition. In particular, front errors are nonmonotone and no Stefan
mesh/time order is claimed. Shrinking width at fixed resolution can eventually
underresolve the phase interval. A zero measured domain difference at these
times is evidence of insensitivity at floating-point resolution, not a
general truncation-error bound for later times or other parameters.

The zero-latent heat experiment uses an insulated 0.1 m slab, 310 K mean,
8 K amplitude and final time 607.927101854 s (one continuum decay time).
Temporal tests use 32 cells and 20/40/80/160 steps. Spatial tests use
10/20/40/80 cells with exact semidiscrete time evolution. These settings were
specified on day 5; only the Stefan settings were fixed on day 1.

## Measurements and tolerance review

See [day-5 results](../validation/enthalpy-day5.md). The existing nonlinear
stopping tolerance is 1e-10 on volume*(rho*c*width+rho*L); the largest measured
step residual in this campaign is below 3e-13. Maximum conservative cumulative
energy defect is below 7e-13 on the whole-slab scale, including unsampled steps.
This checks algebraic closure separately from front accuracy. No tolerance
sweep or guarantee on arbitrary runs is claimed. The frozen method solves its
own linearized equation accurately but can violate true enthalpy conservation.

## Research provenance and limits

This is an original verification implementation, motivated by the read-only
DTU MMME solidification and implicit-conduction topics. No course files,
assignment parameters, teaching code or private data are included. All imported
modules and their provenance manifests remain unchanged. All data here are
synthetic; the results establish neither material calibration nor physical
experimental validation, and make no new-discretization or performance claim.

The similarity construction and its constant-conductivity and prescribed-wall
limits are discussed by A. N. Ceretani, N. N. Salva and D. A. Tarzia,
*An exact solution to a Stefan problem with variable thermal conductivity
and a Robin boundary condition*, arXiv:1706.06984 (2017), sections 2 and 4:
https://arxiv.org/abs/1706.06984. We specialize to constant conductivity and
a prescribed hot wall and derive the melting sign convention explicitly
above; their study treats solidification and a broader boundary model.
The local heat/flux identity tests check the factors and signs independently.
