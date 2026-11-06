A compressor rotor was accelerated quickly through a blade resonance while casing
probes timed every blade tip pass. Each blade's first bending mode was
crossed by an engine-order excitation of unknown order. The rotor is mistuned,
so each blade crosses at its own speed.

The recordings are in `/app/data/`:

- `rig_sheet.md` covers the rotor, the probe layout, the blade numbering, the
  centrifugal stiffening of the mode and the sensor behaviour.
- `opr.csv` holds the once-per-revolution time stamps.
- `probe_1.csv` … `probe_5.csv` hold the raw blade-pass trigger times at each
  probe.

From these recordings, find:

1. the engine order that excites the resonance;
2. for each of the 24 blades, the natural frequency of its first bending
   mode at the speed where that blade crosses resonance, in Hz;
3. for each of the 24 blades, its modal damping ratio, as a fraction of
   critical damping for an equivalent viscous single-degree-of-freedom mode;
4. for each of the 24 blades, the steady-state resonant tip amplitude, in mm
   (zero-to-peak, circumferential). This is the amplitude the blade would
   settle at if the rotor were held at that blade's resonance speed under
   the same excitation, and it is the input to the fatigue assessment.

Write `/app/output/blade_modal_results.json`:

```json
{
  "engine_order": 0,
  "blades": {
    "1":  {"resonance_frequency_Hz": 0.0, "damping_ratio": 0.0, "steady_resonant_amplitude_mm": 0.0},
    "2":  {"resonance_frequency_Hz": 0.0, "damping_ratio": 0.0, "steady_resonant_amplitude_mm": 0.0},
    "...": {},
    "24": {"resonance_frequency_Hz": 0.0, "damping_ratio": 0.0, "steady_resonant_amplitude_mm": 0.0}
  }
}
```

Blade keys are the blade numbers `"1"` to `"24"`, as defined in the rig sheet.

Requirements:

- `engine_order` must be the correct integer.
- Every blade's `resonance_frequency_Hz` must be within 0.25% of its true value.
- Every blade's `damping_ratio` must be within 25% (relative) of its true value.
- The signed relative damping error `(estimate − true) / true`, averaged over
  the 24 blades, must lie between −7% and +7%.
- Every blade's `steady_resonant_amplitude_mm` must be within 15% (relative)
  of its true value.
- All of these must hold together.

Python 3 with numpy, scipy and pandas is available.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
