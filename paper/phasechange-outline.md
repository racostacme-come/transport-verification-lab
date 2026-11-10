# Historical day-5 synthesis outline

Superseded on day 6 by [paper.tex](paper.tex) and [paper.pdf](paper.pdf).
This outline is retained as planning provenance; see
[manuscript evidence](../docs/manuscript-evidence.md) for the current workflow.

Target: replace paper.tex and paper.pdf on November 10 with an original
manuscript of at most five pages including references. Keep the v0.1.0
thermal synthesis at its release tag. No phase-change release until the new
manuscript is compiled, source/data-bound, checked in CI and every page rendered
and visually reviewed.

Proposed title: **Energy conservation and front accuracy across a regularized
melting interval: a reproducible enthalpy verification study**.

1. **Question and scope (about half a page).** How do exact enthalpy increments
   and old-state frozen apparent capacity differ when a step crosses phase
   endpoints? Established methods; original synthesis and verification suite.
   Constant-property 1-D conduction, no physical experimental validation.
2. **Connected methods (about one page).** H(T), inverse, finite-volume
   backward Euler, piecewise-linear convex residual, safeguarded Newton,
   C++17 assembly/LDL kernels and explicit Python interface. Define signed
   cumulative true-energy defect and distinguish the frozen equation residual.
3. **Independent references and protocol (about one page).** Event traversal
   for uniform heating/cooling, algebraic/semidiscrete/continuum cosine
   references, Stefan similarity with the solid initially at Tm. Give fixed
   synthetic parameters, observation window, front definition, and three
   isolated four-level studies plus domain doubling.
4. **Measured results (about one page).** Tables generated from
   results/enthalpy-baseline.json, enthalpy-comparison.json and
   enthalpy-refinement.json. Use heat orders 0.9853/0.9926/0.9963 and
   1.9830/1.9958/1.9989. Finest conservative front errors at 30/60/90/120 s:
   0.066729/0.042863/0.057387/0.213862 percent; frozen:
   1.542954/1.052027/0.830170/0.617133 percent. Retain the case where the
   coarser time step gives a smaller early conservative front error, the
   mesh-1/30 s frozen advantage and the closed-cycle transient energy defects.
5. **Discussion, reproducibility, limitations and references (remaining space).**
   Distinguish implementation agreement, exact discrete evolution, continuum
   error, finite-interval regularization and physical validity. Report domain
   insensitivity only for the observed window. Explain error cancellation;
   no universal superiority or Stefan convergence-order claim. State the
   native backend shares Python safeguarding. Give exact reproduction CLI,
   versioned source hashes, actual hosted checks and release gate status.

Reference candidates already checked: Voller, Cross and Markatos (1987),
DOI 10.1002/nme.1620240119; Ceretani, Salva and Tarzia (2017),
https://arxiv.org/abs/1706.06984, especially sections 2 and 4. Do not attribute
novelty to these established numerical formulations or invent an affiliation.

Next manuscript work must generate tables from source-bound data, add a
stale-table check, update suite.py paper inputs for all claimed experiments,
compile, reject undefined references/layout errors, count all pages including
references, and inspect every locally and independently CI-rendered page.
