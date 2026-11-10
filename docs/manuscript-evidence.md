# Phase-change manuscript evidence

The current synthesis replaces the root thermal paper; the original remains
available unchanged at the v0.1.0 tag. No imported papers or source snapshots
are edited. The question is energy conservation versus front accuracy across
a regularized phase interval, with an intentionally frozen comparator.

## Reproduce and build

After installation and the native build described in README:

```sh
python -m phasechange --output validation/enthalpy-baseline.json
python -m phasechange --study comparison --check-against results/enthalpy-comparison.json --output validation/comparison-python.json
python -m phasechange.refinement --check-against results/enthalpy-refinement.json --output validation/refinement-python.json
python manuscript.py check
python manuscript.py check-report
```

Repeat comparison and refinement with `--backend native` and distinct output
paths. Both OS jobs do this. To intentionally change numerical sources, rerun
and review affected experiments before replacing versioned results;
never merely update hashes to make stale results appear current.

```sh
python manuscript.py generate
python suite.py paper --engine tectonic
python manuscript.py check-report
```

pdfLaTeX is also supported and runs twice. Render for visual review with
`pdftoppm -r 120 -png paper/paper.pdf paper/renders/page` after creating that
ignored directory.

## Meaning and limits of the gates

- Three reports bind uniform baseline, comparison and refinement to numerical
  source hashes. Refinement requires the complete current phasechange/native
  source inventory and all five named acceptance flags.
- Missing, nonfinite or unresolved front errors and incomplete observation
  windows fail. Experiment reproduction and numerical tests in CI check the
  calculations; manuscript checks alone do not independently recompute them.
- Five generated TeX fragments contain uniform, slab, heat and front tables
  and measured macros. The generator never edits numerical data. Checks compare
  generated content without modifying it.
- The build manifest covers all numerical sources, three reports, generated
  fragments, manuscript, generator, runner, workflow, manuscript tests,
  snapshot manifest and compiled PDF.
- The report checker requires the exact input set, every hash and an actual
  page count of 1--5 including references. CI verifies the committed artifact
  before independently rebuilding and checking it. Different engines may
  legitimately produce different PDF hashes.
- Build logs reject overfull boxes and unresolved references. Every final
  local and CI page must also be rendered and visually reviewed.

Mutation tests challenge source drift, omitted/new sources, failed gates,
missing/invalid observations, altered measurements, stale tables, missing
manifest inputs, altered PDF bytes and inconsistent page counts. They do not
establish physical validity or detect every incorrect scientific claim.
Human review still checks definitions, units, references, parameter descriptions
and conclusions against the versioned evidence.
