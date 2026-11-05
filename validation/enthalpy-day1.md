# Phase-change day-1 validation

Scheduled-run date: 5 November 2026, Europe/Paris. The run continued after
midnight on 6 November; Git uses actual commit timestamps. November 3-4 is
the requested break, not missing implementation days.

## Executed locally

- `python -m phasechange`: 30 fixed cases; maximum scaled temperature error
  against independent event traversal 0.0; maximum scaled energy defect
  1.3932639823586679e-14 (gate 1e-12). Temperature scale 2 K; volumetric energy
  scale 110160000 J/m3. Raw samples, input and source hashes are versioned in
  `results/enthalpy-baseline.json`.
- `python -m pytest tests/test_phasechange.py -q`: 48 passed, including
  hand-calculated endpoints, inverse/monotonicity, zero-latent limit,
  one-sided derivatives, reference-offset invariance, crossing/reversal,
  explicit invalid-state rejection and versioned result reproduction.
- `python suite.py test`: 304 passed: 246 imported-component and 58 root tests.
  Exact outputs are recorded in `validation/tests.json`.
- `python suite.py native --build-root ../../.native-transport`: both MSVC
  C++17 CTests passed; these are existing component regressions, not new
  phase-change native kernels. Outputs: `validation/native.json`.
- `python suite.py verify`: all 192 Transport snapshot hashes passed.
- `python suite.py experiment`: both prior connected thermal outputs reproduced
  without a tracked result change.
- Workspace audit: all 598 imported hashes across three suites passed; all
  19 original component repositories were clean; 34 previously indexed lecture
  hashes matched. Five new eligible PDFs were reviewed for topics; private
  exercise/report/code changes were excluded. No original/source files changed.

## Evidence boundaries

This baseline verifies constitutive and uniform heat-input accounting only.
For piecewise constant power the update is exact, so the subdivision sweep is
an implementation/roundoff check, not evidence of temporal convergence.
Spatial phase-change accuracy, nonlinear convergence, front errors, new C++
kernels and frozen-capacity comparison remain future acceptance gates.
No physical or experimental validation is claimed.

The existing three-page v0.1.0 synthesis is unchanged. The new synthesis is
scheduled after connected experiments; no replacement paper is claimed today.
Hosted Linux/Windows and manuscript results must be read from the workflow for
this exact commit; local passes do not imply hosted success. Feature work is
published without merging or releasing the unfinished extension.
