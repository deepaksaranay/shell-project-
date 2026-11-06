# Kepler tasks

| task | bundle | upload zip | generator (not in bundle) |
|---|---|---|---|
| rail-creep-curves (current) | `rail-creep-curves/` | `rail-creep-curves.zip` | `rail-creep-curves/authoring/provenance/` |
| surge-vessel-sizing (failed similarity) | `surge-vessel-sizing/` | `surge-vessel-sizing.zip` | `surge-vessel-sizing/authoring/provenance/` |
| blade-tip-timing (failed similarity) | `blade-tip-timing/` | `blade-tip-timing.zip` | `blade-tip-timing/authoring/provenance/` |
| bearing-impedance-identification (earlier draft) | `bearing-impedance-identification/` | `bearing-impedance-identification.zip` | `generator/` |

- Each generator folder has `generate.py`, which regenerates the data and `tests/expected.json`.
- It also has scripts (`check.py`, `alt_*.py`, `ablate.py`) showing that the
  stated tolerances pass correct methods and fail flawed ones.
