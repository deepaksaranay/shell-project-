# Kepler tasks

| task | bundle | upload zip | generator (not in bundle) |
|---|---|---|---|
| blade-tip-timing (current) | `blade-tip-timing/` | `blade-tip-timing.zip` | `generator_btt/` |
| bearing-impedance-identification (earlier draft) | `bearing-impedance-identification/` | `bearing-impedance-identification.zip` | `generator/` |

- Each generator folder has `generate.py`, which regenerates the data and `tests/expected.json`.
- It also has scripts (`check.py`, `alt_*.py`, `ablate.py`) showing that the
  stated tolerances pass correct methods and fail flawed ones.
