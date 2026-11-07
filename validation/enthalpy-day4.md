# Day 4: matched frozen-capacity experiments

November 8, 2026, early run following completed day 3. Actual timestamps are
retained; this work is not attributed to the later scheduled 23:00 run. The
cycle remains November 5-11, with the requested break excluded. No release.

Original frozen-capacity solver, explicit reuse of the native C++17 kernels,
common CLI, signed cumulative ledger and unambiguous front extraction now
connect the two algorithms on exactly matched inputs. The comparator's
right-sided endpoint convention is documented, including cooling asymmetry.

Local Windows checks:

- 45 new tests cover independent scalar and two-cell equations, endpoint
  crossings, zero-latent ties on both backends, invalid input/state failures,
  hidden transient defects, unresolved fronts and corrupted evidence rejection.
- Full suite: 405 Python tests (246 component + 159 root); three MSVC CTests.
- 30 uniform, five spatial heat/cool and four coarse front comparison cases
  reproduced through Python and native CLI paths. Maximum local sampled
  temperature difference between backends: 1.13687e-13 K.
- All 598 snapshot hashes and 39 indexed lecture hashes match; 19 original
  repositories remain clean. No new eligible course PDFs. Source folders read
  only; no original snapshot, manifest or history changes.
- The baseline source hash is explicitly refreshed for the common CLI change.
  Its numerical JSON is exactly unchanged. Other historical source files and
  result hashes are unchanged. New study sources have a separate full manifest.

## Measured outcomes

The uniform solid/mushy cases retain ties. A single 50 s heating step from
305 K with 2.7e6 W/m3 crosses the entire interval: the frozen method is 100 K
above the independent event reference. The corresponding closed heat/cool
cycle returns to the correct initial temperature with zero final energy
defect, despite large intermediate errors. All-step maxima expose this.

For the spatial heat/cool experiment, the maximum sampled temperature
difference between methods is nonmonotone with refinement:

| dt (s) | Method difference (K) | Frozen maximum cumulative scaled energy defect |
| ---: | ---: | ---: |
| 10 | 69.4535 | 0.980392 |
| 5 | 73.1832 | 0.980392 |
| 2.5 | 74.9719 | 0.980392 |
| 1.25 | 68.6038 | 0.848771 |
| 0.625 | 34.8213 | 0.442114 |

These differences have no independent spatial reference yet. The conservative
maximum cumulative scaled energy defect over these runs is 1.961e-14.

For N=100 and a fixed 1 K melting interval, the 120 s front positions are:

| dt (s) | Enthalpy front (mm) | Frozen front (mm) | Frozen max cumulative scaled defect |
| ---: | ---: | ---: | ---: |
| 0.5 | 8.572968 | 8.468678 | 0.00329943 |
| 0.25 | 8.574811 | 8.537384 | 0.00139610 |
| 0.125 | 8.575734 | 8.552789 | 0.000564576 |
| 0.0625 | 8.576195 | 8.562558 | 0.000456071 |

Every front is resolved at all four observation times. Conservative maximum
cumulative scaled energy defect is 1.490e-13. The frozen maximum *step*
defect slightly worsens between the final two levels (0.000274234 to
0.000289163); that unfavorable outcome is retained. Small algebraic residuals
for the frozen equation do not imply a small true enthalpy error.

Full measured states, signed energy data and source hashes:
[comparison JSON](../results/enthalpy-comparison.json). These are integration
and implementation checks, not continuum or physical validation. The day-1
Stefan parameters and day-5 isolated refinement protocol are unchanged.

The current three-page v0.1.0 synthesis is unchanged supporting prior work.
The phase-change synthesis remains scheduled for November 10 after independent
front/heat validation. Hosted validation state is recorded against the commit
in the persistent portfolio; local success alone does not imply CI success.
