A floating-housing bearing test rig has measured a plain fluid-film journal
bearing at eight operating conditions (four shaft speeds × two static loads).
At each condition the housing was shaken twice with a multisine, first along x
and then along y. Before the oil was turned on, a dry baseline pair of tests
was also recorded.

The raw recordings are in `/app/data/`:

- `rig_sheet.md` describes the rig, the sensors and their calibrations, the
  coordinate and sign conventions, and what each test is.
- `test_log.csv` lists the conditions and their record files.
- `baseline_dry_x.csv`, `baseline_dry_y.csv` and `C01_x.csv` … `C08_y.csv`
  are the time histories.

Identify the oil film's rotordynamic stiffness and damping coefficients at
every condition. Use the model given in the rig sheet:
`f = −(K d + C ḋ + M d̈)`, where `d` is housing motion relative to the shaft.
`K`, `C` and the added mass `M` are frequency-independent 2×2 matrices over
the 20–320 Hz band.

## Output

Write `/app/output/bearing_coefficients.json` with one entry per condition
ID from `test_log.csv`:

```json
{
  "C01": {
    "K_N_per_m":  [[Kxx, Kxy], [Kyx, Kyy]],
    "C_Ns_per_m": [[Cxx, Cxy], [Cyx, Cyy]]
  },
  "C02": { ... }
}
```

Use SI units: N/m for K and N·s/m for C. Row index is the force direction
and column index is the displacement direction, as defined in the rig sheet.

## Requirements

- Each of the 8 conditions must be present with both 2×2 matrices.
- At each condition, every K entry (direct and cross-coupled) must be within
  3% of `max(|Kxx|, |Kyy|)` of that condition's true value.
- At each condition, every C entry must be within 3% of `max(|Cxx|, |Cyy|)`
  of that condition's true value.
- All entries at all conditions must meet these limits at once.

Python 3 with numpy, scipy and pandas is available.

You have 8 hours to complete this task.
