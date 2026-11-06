# Blade-tip-timing test facility: data sheet

This sheet covers how every run from the facility is recorded. Each run
directory holds its own `rig.json` with that rotor's numbers.

## Run directory layout

| file | content |
|---|---|
| `rig.json` | Run constants: `n_blades`; `tip_radius_m`; `southwell_coefficient`; `probe_angles_deg` (nominal, one per probe, in probe order); `probe_angle_tolerance_deg`; `nominal_speed_rpm` (`start`, `end`); `nominal_duration_s`. |
| `opr.csv` | Once-per-revolution trigger times in seconds, one header line (`opr_time_s`). |
| `probe_<k>.csv` | Trigger times in seconds for tip-timing probe k (k = 1 … number of probes), one header line (`toa_s`), in time order. |

All channels share one clock. Timing jitter is about 0.3 µs rms on the
probes and smaller on the OPR.

## Rotors

- Each rotor is a bladed disk (blisk). The disk is very stiff compared with
  the blades, so the blades are mechanically uncoupled.
- The mode of interest is each blade's first bending mode. The probes see its
  tip motion in the direction of rotation.
- Rotation stiffens the blades. The natural frequency follows
  `fn(Ω)² = f0² + S·Ω²`, with `Ω` the shaft speed in rev/s and `S` the
  Southwell coefficient given in `rig.json`. S comes from each blade's
  finite-element model, which has been validated on a spin rig.
- Manufacturing scatter (mistuning) gives every blade a slightly different
  `f0` and damping. Each blade also sits at a slightly different static
  angular position (stagger), and the blades untwist slowly with speed.
- The damping is low, typical of integrally bladed rotors.

## Excitation

Stationary flow distortion from upstream vanes and struts excites the blades
at an integer engine order (EO). Each run crosses one first-bending resonance.
The test records do not document the engine order.

## Runs

The facility runs fast speed sweeps, either accelerations or decelerations,
at roughly 150–200 rpm/s. `nominal_speed_rpm` gives the set-points. The
actual speed is not perfectly smooth, and the OPR record gives the true
shaft speed.

## Sensors

- **OPR pickup:** the pickup on this facility is marginal. In a run it can
  miss a few revolutions and trigger spuriously a few times between
  revolutions. `opr.csv` is the raw trigger list, not a clean
  one-per-revolution record.
- **Tip-timing probes:** optical probes mounted in the casing. Their angles
  are measured in the direction of rotation from probe 1, and each is
  installed to within the stated tolerance.
- **Trigger faults:** the trigger electronics occasionally miss a blade pass
  and now and then fire on noise. A probe's fibre connection can also be
  intermittent for part of a run.

## Blade numbering

At the instant of an OPR pulse, the tip of **blade 1** is aligned with probe
1 (to within its static stagger). Blade numbers increase in the order the
blades pass any given probe. Blade 2 reaches probe 1 one blade pitch after
blade 1, and so on up to blade `n_blades`.
