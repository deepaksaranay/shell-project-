I need the thermal response of our stagnation-point heat shield over the entry trajectory: a charring ablator bonded to an aluminium structure. Please compute it with a proper one-dimensional charring-ablation model: in-depth Arrhenius decomposition, pyrolysis gas flowing out through the char, a receding surface, and a surface energy balance that uses the surface-thermochemistry table with the blowing correction. A heat-of-ablation estimate or a pure conduction model isn't good enough for the bondline sizing.

The inputs are in `/app/data/`:

- `/app/data/case.json`: geometry, initial temperature, the decomposition model and its constants, formation enthalpies, how solid properties are blended between virgin and char, the surface energy balance and blowing correction, interpolation rules, the output times and the definitions of the reported quantities.
- `/app/data/virgin_properties.csv` and `/app/data/char_properties.csv`: specific heat, thermal conductivity and emissivity of the virgin material and the char against temperature.
- `/app/data/pyrolysis_gas.csv`: enthalpy of the pyrolysis gas against temperature.
- `/app/data/surface_thermochemistry.csv`: the B'_c table and the wall gas enthalpy h_w as functions of wall temperature and B'_g.
- `/app/data/trajectory.csv`: the unblown heat-transfer coefficient rhoe_ue_CH0, the recovery enthalpy and the incident radiative flux against time.
- `/app/data/substrate_properties.csv`: the aluminium properties against temperature.

The pyrolysis gas leaves through the heated surface without being stored, in thermal equilibrium with the solid. The aluminium is in perfect thermal contact with the ablator and does not let gas through, and its back face is adiabatic. Run from t = 0 to the last output time.

Write `/app/output/heatshield_response.json` in this format:

```json
{
  "snapshots": [
    {"t_s": 20, "recession_mm": 0.0, "surface_temperature_K": 0.0, "bondline_temperature_K": 0.0,
     "pyrolysis_front_depth_mm": 0.0, "char_front_depth_mm": 0.0, "mass_loss_kg_per_m2": 0.0}
  ],
  "peak_bondline_temperature_K": 0.0
}
```

There is one snapshot for each output time in `case.json`, with the quantities defined there. `peak_bondline_temperature_K` is the highest bondline temperature over the whole run.

Requirements:

- Recession within 3% (plus 0.01 mm) of the converged solution at every output time.
- Surface temperature within 12 K, and bondline temperature within 3 K, at every output time.
- Pyrolysis-front and char-front depths within 3% (plus 0.1 mm).
- Mass loss within 2% (plus 0.01 kg/m²).
- Peak bondline temperature within 3 K.
- The file is a regular JSON file with exactly the keys shown above: top-level `snapshots` and `peak_bondline_temperature_K`, and `t_s`, `recession_mm`, `surface_temperature_K`, `bondline_temperature_K`, `pyrolysis_front_depth_mm`, `char_front_depth_mm` and `mass_loss_kg_per_m2` in each snapshot. Every output time appears once, and every value is a finite number.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
