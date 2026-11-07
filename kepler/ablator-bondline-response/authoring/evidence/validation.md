# Validation record

## Code checks

These checks used modified data sets built in the authoring scratch area.

**Constant-flux slab with an insulated back.** The setup used constant
properties, no decomposition, no recession, zero emissivity and a
negligible substrate. Comparison with the analytic series solution:

| t (s) | surface, exact | generator | solve.py | bondline, exact | generator | solve.py |
|---|---|---|---|---|---|---|
| 40 | 2127.47 | 2127.32 | 2127.39 | 300.130 | 300.136 | 300.133 |
| 100 | 3189.49 | 3189.26 | 3189.37 | 335.990 | 336.049 | 336.021 |

**Recession only, and no decomposition.** The generator and `solve.py`
agree to within 2.3 K at the surface and 0.12 K at the bondline. This
comparison exposed two generator bugs, both fixed before the answer key
was made:

- an upwind index in the grid-motion transport;
- the fibre mass missing from the grid-transported mass.

## Grid convergence of the answer key (generator)

| nodes | peak bondline (K) | recession at 600 s (mm) | mass loss at 600 s (kg/m²) |
|---|---|---|---|
| 120 | 426.27 | 3.511 | 3.991 |
| 240 | 425.75 | 3.511 | 3.986 |
| 480 (key) | 425.49 | 3.511 | 3.983 |

The 240-node run differs from the key by at most 0.30 K (bondline), 0.027 mm
(pyrolysis front) and 0.3% (mass loss).

The convergence is first order: the peak changes by 0.52 K from 120 to 240
nodes and by 0.26 K from 240 to 480. So the key's own discretisation error
is about 0.26 K, which is 9% of the 3 K limit. The independent submission
(440 uniform cells) gives 425.50 K.

## Verifier emulation

The emulation used the same `test_outputs.py` with pytest.

| output | result |
|---|---|
| no file | 7 errors (reward 0) |
| symlink to the answer key | 7 errors (0) |
| copy of the answer key | 7 errors (0) |
| oracle (`solution/solve.py`) | 7 passed (1) |
| independent submission | 7 passed (1) |
| generator, 60 uniform / 100 uniform / 60 graded nodes | 7 passed (1) each |
| no blowing correction | 7 failed (0) |
| pyrolysis-gas heat pick-up ignored | 6 failed (0) |
| no recession | 7 failed (0) |
| B′_c at B′_g = 0 | 1 failed (0) |
| 12-node uniform grid | 6 failed (0) |

## Worst errors, as a share of each limit

| quantity | oracle | independent submission |
|---|---|---|
| recession | −0.0003 mm (0%) | +0.0010 mm (3%) |
| surface temperature | +0.72 K (6%) | −0.61 K (5%) |
| bondline temperature | +0.62 K (21%) | −0.008 K (0%) |
| pyrolysis front | +0.061 mm (17%) | +0.010 mm (5%) |
| char front | +0.002 mm (1%) | +0.016 mm (8%) |
| mass loss | +0.0024 kg/m² (4%) | +0.0013 kg/m² (6%) |
| peak bondline | +0.53 K (18%) | +0.006 K (0%) |

The largest errors of the shortcut models are 2× (B′_c at B′_g = 0) to 30×
(no recession) the limit; `authoring/provenance/margins.py` prints them.
