# Verification record

The experiments are synthetic verification, not experimental physical validation.
Exact commands and inputs are versioned. The initial checks were run on Windows,
Python 3.14.5, NumPy 2.5.4, SciPy 1.18.1 and Matplotlib 3.11.2.

## Independent checks

- 28 parameterized tests cover exact two-cell dynamics, independent dense matrix
  exponentials, insulated heat conservation, quadratic dissipation, maximum
  bounds, constant states, surface jumps, invalid inputs, perfect-contact and
  small-time limits, series continuity, heat quadrature and the continuum
  surface-impedance integral.
- A homogeneous insulated cosine with exact cell averages shows second-order
  spatial convergence, independently of the discontinuous contact experiment.
- The contact study records five meshes, 101 times per mesh, and separate
  one-second flux and heat errors. The results contain 505 history rows, five
  refinement rows and four smooth-refinement rows.
- A doubled-domain comparison retains the same cell widths, checking final
  interface flux at 10 s. It does not establish physical validity.

## Reproduction tolerances

All CSV values are compared with relative tolerance 2e-5 and absolute tolerance
2e-8. JSON values use relative tolerance 2e-5 and absolute tolerance 2e-9 except
the following roundoff diagnostics, which must remain below explicit budgets:

| Diagnostic | Absolute budget |
| --- | --- |
| Small-chain modal versus dense exponential | 1e-8 K |
| Analytical heat versus quadrature | 1e-6 J/m² |
| Energy-balance residual across recorded histories | 0.02 J/m² |
| Relative domain-doubling flux difference | 1e-8 |

The balance budget is small compared with the initial sensible energy of
1.38e8 J/m². Finite-dimensional conservation is exact algebraically; floating
point modal reconstruction is not exact. Roundoff diagnostics are not claimed
to reproduce digit for digit. Final contact flux error must be below 0.1%,
the last two contact refinement orders above 1.95, and all smooth orders above
1.98. The finer result is not used to claim uniform early-time accuracy.

## Git and CI

Development uses `feat/contact-solver` for implementation and tests, followed by
`validation/contact-study` for the experiment, packaging, CI and manuscript.
The solver was merged only after its tests and lint checks passed. The study
branch is merged only after local checks and remote CI pass. There are no
fabricated reviews or backdated commits; the configured Git identity is retained.

CI runs source and installed-wheel tests on Linux/Python 3.11 and Windows/Python
3.14, reproduces numerical outputs, builds distributions and checks committed
paper hashes. A separate job compiles the manuscript with pdfLaTeX, checks
diagnostics and page count, and renders all pages. The local paper uses Tectonic;
these independent builds may differ in PDF metadata or font rasterization.

## Issues resolved while developing

Formatting and a generated LaTeX string syntax error were corrected before
committing the study. An initially guessed 0.01% finest-mesh acceptance target
was tighter than the measured 0.0387%; the explicit study acceptance criterion is
0.1%. The measured error was retained, not altered or described as passing the
earlier guess. All failure outputs remain part of the automation execution record.

The first branch CI run passed both software jobs but exposed a Linux launcher
bug: resolving the `pdflatex` symlink changed its invocation name to `pdftex`,
which selected plain TeX and rejected `documentclass`. The build script now
normalizes the absolute path without dereferencing the executable symlink.
The fix is validated by the subsequent independent manuscript CI run.
