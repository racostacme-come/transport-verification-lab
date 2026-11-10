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

Day 3 adds C++17 residual/Jacobian/LDL kernels and an explicit Python binding.
360 local Python checks and three native CTests pass. The 36 assembly/state
cases and ten-step heat/cool trajectory agree within the declared scales;
see [day-3 evidence](validation/enthalpy-day3.md) and [native API](docs/enthalpy-native.md).
Hosted results are recorded in the persistent portfolio. No merge or release.

Day 4 adds the frozen-capacity comparator, common Python/native CLI, signed
cumulative energy ledger and 30 uniform/five slab/four coarse front experiments.
405 Python tests and three native CTests pass locally; see
[day-4 evidence](validation/enthalpy-day4.md). Ties, nonmonotone refinement and
intermediate defects hidden by final cancellation are retained. The work was
performed early on November 8, with actual commit timestamps preserved.

Day 5 adds independent sharp Stefan and zero-latent heat references, isolated
time/mesh/width/domain studies, original reference identity tests and the
synthesis outline. Both Python and C++ backends reproduce the campaign.
The run began early November 9; actual timestamps are retained. Do not repeat
this milestone at the later scheduled run. See [evidence](validation/enthalpy-day5.md).

Next: November 10 original <=5-page phase-change synthesis with generated,
source-bound tables, compilation and every-page visual review. November 11
release only after every acceptance gate passes. Hosted day-5 status is
recorded in the external portfolio; no merge or release on day 5.
Hosted evidence is recorded in the portfolio. No merge/release before all gates.

## Day 6: synthesis and manuscript evidence gates

November 10 scheduled contribution, continued after midnight with actual Git
timestamps. Original four-page synthesis at [source](paper/paper.tex) and
[PDF](paper/paper.pdf); all tables and measured macros generated from the
source-bound uniform/comparison/refinement reports. Strict source, table,
PDF input-hash and page checks run in CI. Fourteen mutation tests challenge
stale/missing scientific evidence. See [day-6 record](validation/enthalpy-day6.md).

The original v0.1.0 synthesis remains at its tag. Numerical sources, fixed
parameters, original repositories and imported snapshots are unchanged.
Next: November 11 audit all acceptance gates, review the tested branch and
release only if complete; then select a Vibration question from inspected
ATSA/MLSP topics. No merge or release is claimed by the day-6 manuscript work.
