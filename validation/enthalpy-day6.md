# Day 6: original synthesis and evidence gates

Recorded 2026-11-10T23:11:50.244935+00:00. November 10 contribution, continued past midnight;
actual timestamps retained. No merge or release.

## Contribution

Original four-page synthesis in [LaTeX](../paper/paper.tex) and
[PDF](../paper/paper.pdf), including references. All four tables and the
measured macros come from three source-bound JSON reports. The paper connects
uniform endpoint traversal, conservative/frozen finite-volume steps,
Python/C++ assembly, independent heat references and all eleven front settings.
It retains unfavorable and nonmonotone results and makes no physical-validation,
new-method or performance claim. The prior thermal paper remains at v0.1.0.

New manuscript.py rejects stale numerical source hashes, incomplete source
inventories, failed acceptance flags, missing/nonfinite/unresolved front
observations, stale tables, altered PDFs and incomplete build manifests.
The build report binds all 34 required input hashes and the actual
page count. Both OS jobs verify the committed manuscript. The paper job checks
it before independently compiling and checking the new pdfLaTeX artifact.

## Checks actually completed locally

- 445 Python tests passed: 246 component, 199 root (14 new manuscript mutation
  cases). Raw outputs: validation/tests.json.
- Three MSVC C++17 CTests passed via suite.py native; raw build/test output:
  validation/native.json. Numerical kernels were unchanged.
- Uniform baseline reproduced byte-equivalent parsed JSON. Full Python
  comparison (30 uniform, five slab, four coarse fronts) and refinement
  (eleven matched settings, both methods) reproduced against versioned reports.
  All five refinement acceptance flags passed. Reproduction files remain
  locally under ignored build/day6-*.json. Hosted jobs repeat both backends.
- All 34 build inputs verified. Tectonic compile passed: four pages, no
  overfull boxes, undefined references or citations. All four final pages
  rendered with Poppler and visually inspected: legible equations/tables,
  aligned content, complete references, no clipping or orphan final page.
- Initial drafting caught an overflowing equation/command paragraph and an
  orphan reference; these were corrected before the final compile/review.
- Reference metadata checked at the publisher DOI 10.1002/nme.1620240119
  and arXiv:1706.06984. The paper distinguishes the cited broader formulations
  from its independently checked melting specialization.
- All 598 imported hashes across three suites match; all 19 original
  repositories are clean. All 39 indexed lecture hashes match and no new
  eligible lecture PDFs were found. Course sources were read-only. Private
  inventory stays outside the public repository.

## Publication and next step

Local evidence above precedes publication. The exact published commit, hosted
run URL, downloaded artifact digests and independent CI-PDF review are recorded
in the persistent major-projects portfolio and daily evidence directory after
CI completes. This document does not predict those outcomes. GitHub Actions
on this branch exposes the actual status. Day 7 (November 11) reviews all nine
acceptance gates, merges only tested content, and releases only if complete.
Then inspect ATSA/MLSP topics to plan the next Vibration cycle. Preserve this
project if any gate is blocked.
