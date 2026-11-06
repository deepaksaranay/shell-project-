# Validation record

## Verifier emulation

The verifier was emulated with the same test.sh and pytest 8.4.1 +
pytest-json-ctrf 0.3.5. The program under test runs as uid 65534, and the
truth directory is mode 700, root-owned.

| case | outcome | reward |
|---|---|---|
| nop (no /app/btt_analyze.py) | 16 errors | 0 |
| a program that globs and reads `tests/truth/*.json` | finds nothing, crashes; 16 failed | 0 |
| oracle (solution/solve.py) | 16 passed, 4 runs in about 15 s | 1 |

## Reference on the shipped runs

These were run with the pinned numpy 2.3.3, scipy 1.16.2 and pandas 2.3.3.

| run | EO | df max | dz max | dz mean | da max |
|---|---|---|---|---|---|
| example | 6/6 | 0.00028 | 0.113 | +0.005 | 0.065 |
| rotor_a | 7/7 | 0.00028 | 0.085 | +0.019 | 0.059 |
| rotor_b | 5/5 | 0.00029 | 0.058 | −0.003 | 0.040 |
| rotor_c (run-down) | 8/8 | 0.00021 | 0.079 | +0.005 | 0.066 |

## Re-seeded copies

`seed_sweep.sh`, with seed offsets 101, 202 and 303, gives 12 more runs.
The reference passed all 12:

- engine order correct on every run;
- df within 0.00044;
- dz max 0.177;
- dz mean within ±0.014;
- da max 0.096.

## Flawed methods on the shipped runs

**Steady-state SDOF fit plus measured peak.** It fails every run:

| run | df | dz mean | da max |
|---|---|---|---|
| example | 0.0032 | +0.153 | 0.246 |
| rotor_a | 0.0037 | +0.190 | 0.342 |
| rotor_b | 0.0033 | +0.073 | 0.219 |
| rotor_c | 0.0029 | +0.099 | 0.230 |

**No-Southwell variant.** Its mean damping error is:

| run | dz mean | result |
|---|---|---|
| example | +0.145 | fails |
| rotor_a | +0.103 | fails |
| rotor_b | +0.244 | fails |
| rotor_c | +0.055 | passes this run |

All four runs must pass, so the variant fails overall.

**Skipping the OPR repair.** Blade assignment breaks and the fit does not
converge.
