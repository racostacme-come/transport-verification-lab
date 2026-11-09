# Day 5: independent heat and Stefan verification

Recorded 2026-11-09T02:12:34.848934+00:00. Early November 9 run; actual timestamps retained.

## Scope and reproducibility

All three prior consolidations remain released. The active Transport cycle
continues November 5-11, excluding the requested November 3-4 break. No merge
or release today. New synthesis outline: `paper/phasechange-outline.md`;
the unchanged three-page v0.1.0 synthesis is supporting prior work only.

Run `python -m phasechange.refinement` and repeat with `--backend native`.
Use `--check-against results/enthalpy-refinement.json` for source/data checks.
The workflow reproduces Python and native campaigns on Linux and Windows.
Hosted outcome is recorded in the persistent portfolio after publication.

## Fixed front studies

Every entry is the maximum absolute relative front error over 30/60/90/120 s.
Shared finest settings are N=800, dt=0.0625 s, width=0.125 K, length=0.1 m.
Mesh/time/width variations change only the named parameter; domain doubling
retains dx. Full observed fields, signed energy and all inputs are in the JSON.

| Case | N | dt (s) | width (K) | Enthalpy error (%) | Frozen error (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| refined | 800 | 0.0625 | 0.125 | 0.213862 | 1.542954 |
| mesh-0 | 100 | 0.0625 | 0.125 | 2.861279 | 3.068805 |
| mesh-1 | 200 | 0.0625 | 0.125 | 0.839380 | 1.032510 |
| mesh-2 | 400 | 0.0625 | 0.125 | 0.625741 | 1.311259 |
| time-0 | 800 | 0.5 | 0.125 | 0.217875 | 11.078534 |
| time-1 | 800 | 0.25 | 0.125 | 0.214978 | 7.142458 |
| time-2 | 800 | 0.125 | 0.125 | 0.214434 | 2.409037 |
| width-0 | 800 | 0.0625 | 1.0 | 1.583802 | 2.223497 |
| width-1 | 800 | 0.0625 | 0.5 | 0.898116 | 1.511194 |
| width-2 | 800 | 0.0625 | 0.25 | 0.277012 | 1.551998 |
| double-domain | 1600 | 0.0625 | 0.125 | 0.213862 | 1.542954 |

At refined settings both methods meet the 2% front target. Their four
conservative errors are 0.066729%, 0.042863%, 0.057387%, 0.213862%; frozen
errors are 1.542954%, 1.052027%, 0.830170%, 0.617133%. The sharp-interface
reference has lambda=0.3064239053612112. Finite-interval and front-interpolation
errors remain in the measured deviations; no Stefan convergence order claimed.

Nonmonotone/unfavorable outcomes are retained: at 30 s, mesh-1 frozen error
is 0.177702% versus 0.699333% for enthalpy. With N=800 and width=0.125 K,
the conservative 30 s error increases as dt decreases (0.030316% to 0.066729%).
These error cancellations preclude a universal-superiority claim. Doubling
the domain gives zero measured front and common-cell temperature differences
for both methods at the four sampled times; this is no general domain bound.

## Heat limits and algebraic closure

- Temporal orders: 0.98529697, 0.99256946, 0.99626459.
- Spatial orders: 1.98300376, 1.99576739, 1.99894287.
- Temporal error uses exact semidiscrete evolution as its reference; the spatial
  error uses exact semidiscrete evolution versus continuum cell averages.
- The independent incidence-matrix exponential and differential heat/Stefan
  flux identities check the analytical oracle; no solver/material imports.
- Maximum heat algebraic error: 2.842171e-13 K.
- Maximum step residual (both methods): 2.976727e-13 (stopping tolerance 1e-10).
- Maximum enthalpy cumulative true-energy defect: 6.433068e-13 on volume*(rho*c*width+rho*L).
- Maximum frozen cumulative true-energy defect: 4.247508e-02 on volume*(rho*c*width+rho*L).

Residual/ledger closure is distinct from continuum accuracy. No tolerance
sweep, performance claim or physical experimental validation is asserted.

## Local verification and provenance

431 Python tests passed (246 component and 185 root; 26 new day-5 tests).
Three MSVC C++17 CTests passed. Both complete local campaigns passed all
acceptance gates and reproduced within declared tolerances.
The final local suite report is `validation/tests.json`; native build/CTest
output is `validation/native.json`. The day-5 additions comprise 26 new
tests, including extreme similarity scales and rejection of missing/failed
front evidence. All 598 imported file hashes across the three suites pass;
all 19 originals are clean and all 39 indexed lecture hashes match. No new
eligible lecture PDFs. No imported code/manifests or original histories changed.

The comparison manifest explicitly adds the new reference and campaign
sources; all prior comparison numerical values are exactly unchanged. The
day-1 fixed protocol remains byte-for-byte unchanged. Current material inputs
are original synthetic values; course files remain read-only and unpublished.

Next: November 10 generated paper tables, original <=5-page phase-change
synthesis, compile and visual review of every page locally and in CI.
