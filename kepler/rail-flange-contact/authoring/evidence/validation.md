# Validation record

## Verifier runs

`authoring/run_verifier_checks.sh` copies `tests/`, points the result path
at a scratch file, and runs the unchanged `test_outputs.py` once per file
in `verifier_inputs/`. Every file there is in exactly the output schema:

- `oracle.json`: the output of `solution/solve.py`;
- `independent_submission.json`: the output of
  `authoring/provenance/independent_submission.py`, which writes the schema
  itself;
- `generator_el0.4.json` and `baseline_*.json`: generator and shortcut
  results exported by `authoring/provenance/export_for_verifier.py`.

The full pytest output is in `verifier_runs.log`.

| output | result |
|---|---|
| no file | 6 errors (reward 0) |
| symlink to the answer key | 6 errors (0) |
| copy of the answer key | 6 errors (0) |
| oracle (`solution/solve.py`) | 6 passed (1) |
| independent submission | 6 passed (1) |
| generator at 0.4 mm | 6 passed (1) |
| creep forces left out of the equilibrium | 6 failed (0) |
| horizontal patches | 5 failed (0) |
| one rolling radius for both patches | 6 failed (0) |
| 1 mm grid | 2 failed (0) |

## Per-case errors against the key

Each line gives the case, then for each patch the contact-angle error (deg), the normal-force error (N) and the creep-force vector error (N).

**Oracle:**

```
A  shift +0.0000  tread: ang -0.001 N +13 F 7  flange: ang -0.032 N -12 F 2
B  shift -0.0004  tread: ang +0.004 N -149 F 50  flange: ang -0.005 N +111 F 139
C  shift -0.0002  tread: ang +0.004 N -165 F 65  flange: ang +0.004 N +129 F 126
D  shift +0.0002  tread: ang +0.001 N +37 F 30  flange: ang -0.088 N -38 F 9
E  shift -0.0003  tread: ang +0.002 N -133 F 64  flange: ang -0.017 N +102 F 125
F  shift -0.0003  tread: ang +0.004 N -116 F 16  flange: ang -0.023 N +65 F 132
worst share of limit: {'shift': '4%', 'angle': '29%', 'N': '21%', 'F': '22%'}
```

**Independent submission:**

```
A  shift +0.0001  tread: ang +0.003 N +14 F 3  flange: ang -0.042 N -7 F 3
B  shift +0.0001  tread: ang -0.003 N +54 F 9  flange: ang -0.037 N -10 F 52
C  shift +0.0001  tread: ang +0.001 N +32 F 23  flange: ang -0.037 N -25 F 107
D  shift -0.0002  tread: ang +0.001 N -25 F 80  flange: ang -0.045 N +85 F 28
E  shift +0.0000  tread: ang -0.006 N +47 F 12  flange: ang -0.037 N -10 F 41
F  shift +0.0001  tread: ang +0.004 N +66 F 50  flange: ang -0.037 N -63 F 64
worst share of limit: {'shift': '2%', 'angle': '15%', 'N': '11%', 'F': '13%'}
```

## Code checks

The contact-search definition was changed from local-minimum basins to fixed wheel regions, because linear-interpolation noise in the nearly conformal tread region created spurious extrema. All three codes use the region definition.
