# Throat-erosion ballistics

## Task

One static firing of a small BATES-grain solid rocket motor is given as a 2 kHz record of gauge chamber pressure and
un-zeroed load-cell thrust, plus motor, propellant and nozzle data. The graphite throat erodes during the burn and only its pre- and
post-fire diameters are known. The solver must reconstruct the delivered characteristic velocity $ c^* $, the nozzle jet efficiency, the
burning-rate law $ r = r_4 (p_c/4\,\mathrm{MPa})^n $, the total impulse and the throat-diameter history.

The intended route combines two inverse problems:

1. **Thrust-coefficient inversion.** $ F/p_c $ depends on the throat area both directly and through the area ratio
   $ A_e/A_t $ (isentropic $ p_e/p_c $). Inverting it sample by sample gives $ A_t(t) $ up to the unknown jet efficiency. The efficiency
   is calibrated in the early burn, before the throat insert has heated enough to erode, where $ A_t $ equals the measured pre-fire value.
   The post-fire measurement then serves as a check.
2. **Single-firing burn-rate reduction.** With $ A_t(t) $, the overall mass balance gives $ c^* $. The discharged mass, corrected for
   gas stored in the growing free volume, gives the web regression $ w(t) $ through the BATES geometry. Its derivative against
   $ p_c $ gives $ n $ and $ r_4 $.

## Difficulty

Propulsion test engineers at solid-motor manufacturers, propellant-development labs and university or amateur
rocketry test stands do this reduction after every development firing. The classic single-firing burn-rate method (Nakka;
mass-balance methods in the DRDO literature) assumes a constant throat. With a throat that opens up by about 7 % in
diameter, that assumption collapses the pressure range and gives a meaningless exponent. The usual practitioner fix is to
interpolate the throat linearly in time between the pre- and post-fire measurements, or in proportion to the cumulative pressure. Both
ignore the thermal onset delay of throat erosion. They bias the throat history by about 0.2 mm and $ n $ by about 0.04,
and the error does not show in the pressure or thrust fits. Recovering the throat history from the thrust coefficient is the
subject of recent work (time-resolved throat reconstruction from coupled thrust and pressure, demonstrated so far on hybrids). The
solver also has to handle sensor details a test engineer would: gauge versus absolute pressure, the load-cell zero offset, igniter gas
in the mass balance, and subsonic nozzle flow during start-up and tail-off.

The data are synthetic. They come from a lumped 0-D internal-ballistics simulation with chamber gas storage, a choked/subsonic nozzle flow
function, pressure-driven throat erosion after a thermal onset delay, divergence and jet-efficiency thrust losses, a load-cell
offset and Gaussian sensor noise (`authoring/provenance/generate_data.py`). Motor size, propellant properties (AP-composite-like
$ \gamma $, density, $ c^* $, $ n $), erosion rate and noise levels are in realistic ranges for a 75–100 mm research motor. Real
motors add effects that are not modelled here: erosive burning, 2-D burnback, nozzle heat loss and transducer dynamics.

Constant answers cannot pass: every graded quantity depends on the specific record.

## Verification

Evidence is in `authoring/evidence/`. The reference solution (`solution/solve.sh`) scores reward 1
(`oracle_test_log.txt`) and the do-nothing run scores 0 with a valid CTRF report (`nop_test_log.txt`).
The tolerances were calibrated against an independent implementation, `independent_check.py`. It fits a forward 0-D quasi-steady
simulation with a free piecewise-linear throat history to the pressure and thrust records by least squares, with no
thrust-coefficient inversion. It passes all tests (`independent_test_log.txt`; it does not compute total impulse, so the reference value was
used for that key).

| Quantity | Truth | Reference solution | Independent forward fit | Tolerance |
|---|---|---|---|---|
| `c_star_m_per_s` (m/s) | 1492.0000 | 1496.4608 | 1476.0989 | 1.5 % |
| `jet_efficiency` (-) | 0.9520 | 0.9520 | 0.9498 | 0.01 |
| `burn_rate_exponent` (-) | 0.3600 | 0.3618 | 0.3613 | 0.02 |
| `burn_rate_at_4MPa_mm_per_s` (mm/s) | 6.4585 | 6.4404 | 6.5137 | 2 % |
| `total_impulse_N_s` (N·s) | 7134.5101 | 7134.7729 | n/a | 1 % |
| throat diameter at 1.0 s (mm) | 18.0482 | 18.0486 | 18.0692 | 0.05 mm |
| throat diameter at 1.5 s (mm) | 18.2150 | 18.2145 | 18.2369 | 0.05 mm |
| throat diameter at 2.0 s (mm) | 18.3899 | 18.3893 | 18.4143 | 0.05 mm |
| throat diameter at 2.5 s (mm) | 18.5710 | 18.5674 | 18.5933 | 0.05 mm |
| throat diameter at 3.0 s (mm) | 18.7567 | 18.7577 | 18.7824 | 0.05 mm |

Spread of correct-method variants and failure of the common shortcuts (`shortcut_and_variant_results.txt`):

```
# Correct-method variants (should pass)
reference :: 12 passed in 0.02s
no gas-storage term in web regression :: 12 passed in 0.02s
choked flow assumed at all pressures :: 12 passed in 0.02s
igniter gas mass ignored :: 12 passed in 0.02s
heavier throat smoothing (0.6 s) :: 12 passed in 0.02s
steady-window threshold 60 % :: 12 passed in 0.02s
# Shortcuts / errors (should fail)
constant pre-fire throat (classic single-firing method) :: 7 failed, 5 passed in 0.06s
throat linearly interpolated pre->post in time :: 8 failed, 4 passed in 0.04s
throat erosion weighted by cumulative pressure (no onset delay) :: 7 failed, 5 passed in 0.04s
gauge pressure used as absolute :: 6 failed, 6 passed in 0.04s
load-cell zero offset not removed :: 7 failed, 5 passed in 0.05s
divergence factor omitted :: 1 failed, 11 passed in 0.03s
```

The tolerances sit at roughly 1.5–2× the largest deviation among correct methods (independent fit: $ c^* $ −1.07 %; choked-flow
simplification: $ n $ +0.009, $ r_4 $ −0.9 %). The shortcut methods miss the throat history by about 0.2 mm or more, and $ n $ by
0.04 or more.

## Layout

- `instruction.md`, `task.toml`
- `environment/` – image and data (`/app/data/static_fire.csv`, `/app/data/motor.json`)
- `solution/solve.sh` – reference reconstruction
- `tests/` – verifier (`test.sh`, `test_outputs.py`, `expected.json`, `Dockerfile`)
- `authoring/hint.md`, `authoring/provenance/` (data generator), `authoring/evidence/` (logs, results, independent check)
