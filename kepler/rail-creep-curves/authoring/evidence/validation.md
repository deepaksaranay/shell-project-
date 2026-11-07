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

## Independent implementation through the verifier

`authoring/provenance/independent_check.py` solves the steady-state
equations directly with an active-set method, on a 0.45 mm grid shifted by
a quarter element. It scored 10/10 passed (reward 1): area +1.51%, peak
pressure +0.70%. Its force errors, as a share of μN:

| A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|
| 0.47% | 1.07% | 0.79% | 1.01% | 1.06% | 1.75% | 0.86% | 1.14% |

Cases B, F, G and H stopped at the 3000-iteration cap. They were in a
small 2-cycle at one stick/slip boundary element, and their forces are
within the limits shown above.

## Structurally independent submission through the verifier

`authoring/provenance/independent_submission.py` was written separately
from the generator and from `solution/solve.py`. It uses:

- non-square 0.30 × 0.40 mm elements on a grid offset by 0.07 mm;
- zero-padded FFT convolution of the closed-form rectangle kernels;
- a Polonsky–Keer conjugate-gradient normal solver;
- ADMM for each rolling step of the transient scheme.

Its output (`independent_submission_output.json`) scored 10/10 passed
(reward 1): area +1.14% (208.4 mm²), peak pressure +1.87% (830.8 MPa).
Its force errors, as a share of μN:

| A | B | C | D | E | F | G | H |
|---|---|---|---|---|---|---|---|
| 0.29% | 0.69% | 0.50% | 0.66% | 0.67% | 1.35% | 0.72% | 0.77% |

## Anti-cheat emulation

| output | reward |
|---|---|
| no file | 0 |
| symlink to the answer key | 0 |
| copy of the answer key | 0 |
| independent submission | 1 |
| oracle | 1 |
