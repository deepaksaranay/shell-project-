# Provenance

- `environment/data/system.json`: hand-built design case. `tests/system.json` is an identical copy for the verifier.
- `tests/hydraulics.py`: the checking MOC model. It uses integer reaches, wave-speed adjustment and Newton iteration at the pump node.
- `solution/solve.py`: an independent model with space-line interpolation, scipy `fsolve` and Heun speed integration. Comparison on five designs (dt = 0.01 s reference, dt = 0.005 s checking):

| design | reference min head | checking min head | difference |
|---|---|---|---|
| 17.9 m³, 0.445, 0.12 m | 3.899 m | 3.893 m | 0.006 m |
| 18.0 m³, 0.42, 0.10 m | 3.045 m | 3.039 m | 0.006 m |
| 22 m³, 0.40, 0.06 m | 5.969 m | 5.964 m | 0.005 m |
| 24 m³, 0.45, 0.12 m | 9.644 m | 9.638 m | 0.006 m |
| 18.5 m³, 0.45, 0.25 m | 4.699 m | 4.693 m | 0.006 m |

  Max gauge head agreed to within 0.01 m and the gas fraction to within 0.0001 on every design.
- `optimize.py` (coarse scan) and `refine.py` (fine scan) locate the minimum feasible volume, about 16.9 m³ at about 45% initial gas. Run them from the bundle root: `python3 authoring/provenance/refine.py 0.01 16.5,17.0`.
- Time-step convergence of the metrics, for (V, gas fraction, d):

| design | dt (s) | min abs head (m) | gas fraction |
|---|---|---|---|
| 22 m³, 0.40, 0.06 m | 0.02 | 5.984 | 0.7644 |
| | 0.01 | 5.969 | 0.7645 |
| | 0.005 | 5.964 | 0.7645 |
| 30 m³, 0.50, 0.15 m | 0.02 | 13.680 | 0.8459 |
| | 0.01 | 13.639 | 0.8460 |
| | 0.005 | 13.628 | 0.8460 |
| 24 m³, 0.45, 0.12 m | 0.02 | 9.659 | 0.8174 |
| | 0.01 | 9.644 | 0.8175 |
| | 0.005 | 9.638 | 0.8175 |

- Feasible gas window at V = 18.0 m³ (d = 0.10 m, dt = 0.01 s): 0.42–0.46. At 18.5 m³ it widens to 0.41–0.47.
- Inflow orifice: the return upsurge exceeds 105 m when d is 0.20 m or more (18.5 m³ vessel, 0.45 gas fraction).
