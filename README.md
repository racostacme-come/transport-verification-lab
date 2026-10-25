# Transport Verification Lab

Conservative transport, incompressibility, boundary layers, reduced thermal models, startup damping and contact resistance.

This major project integrates 6 previously independent research implementations into one reproducible laboratory. Full Python and C++ source, tests, datasets, documentation, licenses and component papers are included under `modules/`; no submodule downloads are required.

## Start here

Python 3.11+ (CI uses 3.12), Git, internet for initial dependency installation, and a C++17 compiler are required. On Windows use Visual Studio 2022 with Desktop development with C++ and a Windows SDK. On Linux use GCC or Clang. Use a fresh virtual environment.

```sh
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
python suite.py verify
python suite.py install
python suite.py test
python suite.py native
python suite.py experiment
```

The runner stops on install failure and reports every module's test outcome. Tests execute in separate processes to avoid duplicate test-file/import names. `validation/tests.json` contains the actual outputs; `results/integration.json` records the connected experiment. Rerunning overwrites only these generated reports. For Windows builds, if compiler discovery fails, use a Visual Studio developer shell or set `CMAKE_GENERATOR=Visual Studio 17 2022`.

## Included components

| Component | Implementation and original paper |
| --- | --- |
| fluxledger | [Source and usage](modules/fluxledger/README.md), [paper](modules/fluxledger/paper/paper.pdf) |
| staggerlab | [Source and usage](modules/staggerlab/README.md), [paper](modules/staggerlab/paper/paper.pdf) |
| layerlens | [Source and usage](modules/layerlens/README.md), [paper](modules/layerlens/paper/paper.pdf) |
| podguard | [Source and usage](modules/podguard/README.md), [paper](modules/podguard/paper/paper.pdf) |
| heatstart | [Source and usage](modules/heatstart/README.md), [paper](modules/heatstart/paper/paper.pdf) |
| contactscale | [Source and usage](modules/contactscale/README.md), [paper](modules/contactscale/paper/paper.pdf) |

## What is integrated

ContactScale and HeatStart share the same heterogeneous insulated finite-volume chain at zero contact resistance. The experiment checks their assembled coefficients and independently computed matrix-exponential references, then measures backward Euler, Crank-Nicolson and Rannacher errors on that shared grid. A second experiment uses exact cell-average cosine modes to separate spatial continuum error from time-integration error without an eigensolver oracle. Two-cell analytic exchange, nonuniform-grid heat/energy/offset checks, and a sharp forward-Euler positivity example support the handoff. The fluid, boundary-layer and reduced-model modules remain separately selectable through their existing APIs.

This is a modular integration, not a universal multiphysics solver. Component APIs retain their explicit assumptions and units. Not every module is physically coupled to every other module.

## Paper and evidence

Read the new [synthesis paper](paper/paper.pdf) and [LaTeX source](paper/paper.tex).
The component papers are retained as detailed supporting studies. New numerical claims are based on the connected experiment, not on simply adding together earlier performance claims.

Compile with Tectonic or pdfLaTeX:
```sh
python suite.py paper --engine tectonic
```
The command rejects unresolved references, overfull boxes and more than five pages. CI independently builds the paper. [Architecture](docs/architecture.md), [source provenance](suite.json), and [development roadmap](ROADMAP.md) describe scope and future work.

## Provenance and development

Imported snapshots are pinned to original Git commits in `suite.json`, with SHA-256 for every imported file. `verify` detects unrecorded changes. Original repositories and histories remain available; this repository starts a new integration history. Intentional module changes require a documented patch record and manifest refresh. MIT notices are retained for every component; the integration layer is also MIT licensed.

CI runs all component Python tests, ten root numerical checks and two native C++17 CTests on Linux and Windows. `python suite.py native --build-root ../../.native-transport` uses a shorter local path if Windows MSBuild encounters path limits. Native sanitizers and full historical reproduction campaigns are not claimed. Recorded results are synthetic numerical verification, not experimental physical validation.

The 25 October local run passed 246 component tests, ten root tests and both native checks. At 80 cells, the homogeneous continuum mass-norm error is 4.847e-5 with observed order 1.9997. Smooth-mode temporal orders from 32 to 64 steps are 0.9921 (BE), 2.0001 (CN) and 2.0002 (Rannacher). See [validation scope](docs/validation.md), [independent results](results/independent.json), and [provenance repair](docs/provenance-repair.json). Hosted CI status must be read from the corresponding workflow run; local success is not a claim of external success.

## Shared API and reproducible outputs

`integration.experiment()` returns the heterogeneous handoff report. `verification.experiment()` returns independent continuum and temporal refinement rows; `verification.cosine_case(N)` supplies the two models and analytic mode/rate. `python suite.py experiment` writes `results/integration.json`, `results/independent.json` and their LaTeX tables. All material properties and times are explicit in these functions and in the paper. Use `python -m pytest tests/test_independent.py -q` for just the independent checks.

Current DTU MMME heat-transfer topics (explicit positivity and interface-cell matching) motivate these original verification experiments. Private lectures, assignments and solutions are excluded. Course-source folders remain read-only.
