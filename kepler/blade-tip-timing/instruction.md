Our blade-tip-timing facility needs an analysis program for its blisk resonance tests. In each test a rotor is swept quickly through a first-bending resonance while casing probes time every blade-tip pass. The engine order that drives the resonance is not documented, and the rotors are mistuned, so each blade crosses resonance at its own speed.

`/app/data/rig_sheet.md` describes how runs are recorded: the run directory layout, the conventions, the blade numbering, the centrifugal stiffening law and the known sensor faults. `/app/data/example/` is one complete run in that format.

Write the program as a single Python file, `/app/btt_analyze.py`. It is invoked as

    python3 /app/btt_analyze.py <run_dir> <output_json>

It must read one run directory and write a JSON file with this content:

```json
{
  "engine_order": 6,
  "blades": {
    "1": {"resonance_frequency_Hz": 0.0, "damping_ratio": 0.0, "steady_resonant_amplitude_mm": 0.0},
    "2": {"resonance_frequency_Hz": 0.0, "damping_ratio": 0.0, "steady_resonant_amplitude_mm": 0.0}
  }
}
```

There is one entry per blade, keyed `"1"` to `"<n_blades>"`. The fields are:

- `engine_order`: the integer engine order of the resonance crossed in the run.
- `resonance_frequency_Hz`: the natural frequency of the blade's first bending mode at the speed where that blade crosses resonance.
- `damping_ratio`: the blade's modal damping as a fraction of critical, for an equivalent viscous single-degree-of-freedom mode.
- `steady_resonant_amplitude_mm`: the zero-to-peak circumferential tip amplitude the blade would settle at if the rotor were held at that blade's resonance speed under the same excitation. This is the input to the fatigue assessment.

The program will also be run on other runs from the facility in the same format. Those rotors differ in blade count, tip radius, Southwell coefficient, probe count and layout, engine order and speed range, and a run can be an acceleration or a deceleration. It is run with Python 3.13 and only the standard library, numpy 2.3.3, scipy 1.16.2 and pandas 2.3.3. It runs as an unprivileged user with the run directory and output path as its only inputs. It must finish each run within 10 minutes on 2 CPU cores.

Requirements, on the example run and on every other run:

- `engine_order` is the correct integer.
- Every blade's `resonance_frequency_Hz` is within 0.25% of its true value.
- Every blade's `damping_ratio` is within 25% (relative) of its true value, and the signed relative damping error `(estimate − true) / true` averaged over the run's blades lies between −7% and +7%.
- Every blade's `steady_resonant_amplitude_mm` is within 15% (relative) of its true value.

All of these must hold together.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
