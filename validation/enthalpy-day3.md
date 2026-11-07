# Day 3: original native enthalpy kernels

November 7, 2026, following the completed November 6 Python slab contribution.
The November 3-4 break remains excluded; target November 11, repair window
November 12-13. This feature is not a completed phase-change release.

Executed locally on Windows with MSVC 19.44:

- `python suite.py native --build-root ../../.native-transport`: three CTests
  passed, including the new phase-change energy/Jacobian/LDL/failure invariants.
  Checks use explicit failures and remain enabled in Release mode.
- `python -m pytest tests/test_phasechange_native.py tests/test_enthalpy_slab.py -q`:
  56 passed (31 new native-interface checks and 25 existing slab checks).
- `python suite.py test`: 360 passed (246 imported-component, 114 root).
- `python -m phasechange.native_checks`: all 36 assembly/state cases and the
  ten-step connected heat/cool trajectory passed. Maximum scaled differences:
  residual 5.614164825758816e-13; Jacobian 1.4302390011274155e-16;
  state 2.842170943040401e-14. Maximum native scaled energy defect:
  1.0526559048297779e-13. Scales and thresholds: `docs/enthalpy-native.md`.
- Earlier slab numerical JSON values were compared after JSON normalization
  and are exactly unchanged. Its source hash was explicitly refreshed for the
  backend hook. An initial tuple-versus-list audit assertion prevented this
  refresh and the stale-hash test failed; after correcting the audit's
  representation comparison, the targeted and full regressions passed.
- All 598 imported files match their manifests; all 19 original repositories
  are clean; all 39 indexed lecture hashes match. No new eligible lecture PDFs
  found. Course sources, imported snapshots and original histories unchanged.

The explicit native path shares Python safeguarding and ledger logic; kernel
agreement is not an independent nonlinear-algorithm validation. Front accuracy,
discretization studies, frozen-capacity comparison and the new synthesis are
future gates. The current three-page v0.1.0 paper is supporting prior work.
Hosted results are recorded against the feature commit in the portfolio and
workflow artifacts. Local results do not imply external success.
