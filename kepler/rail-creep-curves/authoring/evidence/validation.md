# Validation record

## Solver checks against Kalker's table (circular contact)

| ν | quantity | computed | Kalker table |
|---|---|---|---|
| 0 | C11 = C22 | 3.45 | 3.40 |
| 0.25 | C11 | 4.17 | 4.12 |
| 0.25 | C22 | 3.73 | 3.67 |

At large creepage the force saturates at μN. The naive steady-state block
Gauss–Seidel diverges (spectral radius about 20–100), so the transient-scheme
convex formulation is used instead.

## Grid convergence

Values from `problem.py`, with the reference Richardson-extrapolated from the
two grids:

| element size | contact area (mm²) | peak pressure (MPa) |
|---|---|---|
| 0.5 mm | 209.0 | 827.9 |
| 0.35 mm | 208.1 | 824.2 |
| extrapolated reference | 206.1 | 815.6 |

## Force error as a share of μN

| method | worst force error | other |
|---|---|---|
| exact, 0.35 mm grid | 1.1% | — |
| exact, 0.5 mm grid | 1.6% | — |
| exact, 1.0 mm grid | 3.7% | area +4.8%, pressure +4.4% |
| Hertz ellipse (local curvature) | 3.6% | area +11.6%, pressure −10% |
| linear theory, capped at μN | up to 45% | — |

The limits are 2.5% of μN on forces, 4% on area and 5% on peak pressure.

## Verifier emulation

The emulation used the same test.sh with pytest 8.4.1 and
pytest-json-ctrf 0.3.5.

| case | outcome | reward |
|---|---|---|
| nop | 10 errors | 0 |
| oracle (`solution/solve.py`, 0.4 mm grid, about 11 min) | 10 passed | 1 |

The oracle's errors were area +1.70% and peak pressure +1.58%. Its force
errors, as a share of μN:

| A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|
| 0.48% | 0.90% | 0.71% | 0.85% | 0.90% | 1.81% | 0.93% | 1.13% |
