# Conservative slab verification

Run the November 6 cycle-day-2 experiments from the repository root with
`python -m phasechange.slab_checks`. The versioned [slab history](../results/enthalpy-slab.json)
contains all cell temperatures, enthalpies, liquid fractions, iteration counts,
residuals, signed heat contributions in joules and source hashes.

A 12 mm slab of area 0.003 m² starts at 305 K on 24 equal cells. Its right wall
is insulated. For five 10 s steps the left wall is at 335 K and the uniform
source supplies 2.7 MW/m³. Five more steps reverse the source and set the left
wall to 295 K. The final heating step melts every cell. Cooling returns some
cells to solid. Over those ten steps the maximum scaled cell residual was
8.2841e-14 and the independent maximum scaled energy defect was 1.8061e-15.

With latent heat set to zero, a separate insulated 12 s step damps a single
cosine mode. Relative to the exact backward-Euler amplification for this
discrete operator, the maximum temperature error was 5.6844e-14 K. An initial
internal stopping threshold of 1e-12 stalled at a scaled residual of 2.586e-12
in this case. The solver uses the published 1e-10 cell gate. The independent
1e-9 energy gate also applies to every step.

The 25 new tests cover both endpoint directions, one-cell geometry, a steady
linear profile with balanced wall heat, spatially varying manufactured source,
input preservation on iteration failure, discrete-mode decay and history
reproduction. Local validation passed all 329 Python tests (246 component,
83 root) and the two existing native CTests. The new solver is Python/SciPy.
These results verify the discrete implementation. Continuum refinement, front
accuracy, native phase-change kernels and the new synthesis have separate
milestones. All parameters and results are synthetic.
