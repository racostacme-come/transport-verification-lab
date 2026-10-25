# Validation evidence, 25 October 2026

- Installation rebuilt all six packages, including FluxLedger and StaggerLab C++ bindings, with Python 3.12.2 on Windows.
- `python suite.py test`: 246 component tests (42, 41, 42, 51, 42, 28) and ten root numerical checks passed. Full output: `validation/tests.json`.
- `python suite.py native --build-root ../../.native-transport`: both C++17 invariant executables passed, compiled using MSVC 19.44. Full commands/output: `validation/native.json`. Checks remain active in Release builds (explicit failure branches, not disabled assertions).
- `python suite.py experiment`: measured heterogeneous handoff and independent homogeneous cosine convergence. Generated JSON and LaTeX tables live in `results/`.
- Synthesis compiled with Tectonic; three pages including references; every page rendered and visually inspected. No overfull boxes or undefined references. Build input/output hashes: `paper/build-report.json`.
- All 192 imports were compared byte-for-byte against the pinned original Git blobs. Twenty manifest hashes had described CRLF archive conversion, not those canonical blobs; old/new hashes are recorded in `docs/provenance-repair.json`. The only intentional imported edit pins StaggerLab's default text checkout to LF. Numerical source is unchanged.

- Fresh clones with `core.autocrlf=false` and `true` both passed all 192 byte-integrity checks (`validation/checkouts.json`).
- All 19 original component repositories remained clean (`validation/original-repositories.json`).
- All nine synthesis build-report input/output hashes matched the committed files. Final local PDF pages were visually inspected after the title/layout revision.

## What each check establishes

Implementation: separate assembly agreement, closed-form unequal two-cell exchange, analytic cosine discrete eigenmode, invariant heat/energy/offset tests, and native transport/operator identities.

Discrete error: time stepping versus the analytic discrete solution on a fixed grid. BE approaches first order; CN and fixed two-interval Rannacher startup approach second order for the smooth mode. The discontinuous heterogeneous initial field has its own measured errors and does not establish general asymptotic rates.

Continuum error: exact-in-time finite volumes versus exact cell averages of a smooth homogeneous insulated heat solution. Observed order approaches two. This does not establish an interface-error rate for discontinuous materials or nonzero contact resistance.

Physical validity: no laboratory experiment or material calibration. Positive-coefficient forward Euler is checked locally; negative coefficients are not equated with spectral instability. No sanitizer or full historical reproduction campaign is claimed.

## Publication gates

Linux and Windows CI repeat numerical/native checks; a separate pdfLaTeX job builds and checks the manuscript and uploads it. Release requires successful hosted checks and review of the independently compiled paper. Actual run URLs/status are recorded in the portfolio outside this repository.
