# Transport enthalpy-conservation cycle

Active study: When phase-change steps cross a melting interval, how do conservative enthalpy increments and a frozen apparent-heat-capacity approximation differ in energy balance, temperature and front accuracy?

The requested November 3-4 break is excluded. Day 1 began in the November 5
scheduled run, which continued past midnight. Git records actual timestamps.
Target November 11; earliest finish November 9; November 12-13 only justified
validation/polish. At day 9 retain this project and report any blocker/recovery.

All three original suite consolidations are released. Transport v0.1.0 is
preserved. New study source and acceptance contract: [enthalpy conservation](docs/enthalpy-conservation.md).

| Day | Date | Planned contribution |
| --- | --- | --- |
| 1 | 2026-11-05 | Confirm the question, fixed synthetic benchmark parameters, scaling and acceptance criteria; document architecture; implement constitutive enthalpy/inverse and independent zero-dimensional heating baseline; commit and push. |
| 2 | 2026-11-06 | Implement conservative finite-volume backward Euler with safeguarded nonlinear solve, energy ledger, phase-endpoint failure cases and zero-latent heat reference checks. |
| 3 | 2026-11-07 | Implement original C++17 residual/Jacobian and tridiagonal kernels, explicit Python interface, independent assembly/derivative checks and standalone native invariants. |
| 4 | 2026-11-08 | Integrate the frozen-capacity comparator and common CLI; run original heating/cooling and phase-front experiments with source-bound outputs and parameter sweeps. |
| 5 | 2026-11-09 | Complete independent Stefan and heat references, separate time/mesh/regularization studies, review nonlinear tolerances and reproduce Linux/Windows CI artifacts; begin synthesis outline. |
| 6 | 2026-11-10 | Reproduce selected results, draft the original <=5-page synthesis with generated tables and verified references, compile in CI and locally, render and inspect every page. |
| 7 | 2026-11-11 | Review every acceptance gate, merge only tested content, publish release if complete and select the next coherent Vibration question from inspected ATSA/MLSP topics. |
| 8 | 2026-11-12 | Use only for explicitly justified numerical, portability or manuscript repairs; retain this project. |
| 9 | 2026-11-13 | Finish justified repairs or report precise blocker and recovery plan; continue the same project until acceptance rather than falsely declare completion. |

Day 2 delivered the Python conservative finite-volume step and independent
energy ledger. The [day-2 record](validation/enthalpy-day2.md) documents the
numerical checks. Next, implement the C++17 residual, Jacobian and tridiagonal
kernels and compare them with the Python solver.
