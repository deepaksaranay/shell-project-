# Compressor rotor blade-tip-timing run: data sheet

## Rotor

- Bladed disk with **24 blades**. The tip radius is **247.5 mm**.
- The disk is very stiff compared with the blades, so the blades are
  mechanically uncoupled.
- The mode of interest is each blade's first bending mode. The probes see its
  tip motion in the direction of rotation.
- Rotation stiffens the blades. The natural frequency of this mode follows
  `fn(Ω)² = f0² + S·Ω²`, with `Ω` the shaft speed in rev/s and **Southwell
  coefficient S = 4.2**. S comes from the blade finite-element model, which
  has been validated on a spin rig.
- Manufacturing scatter (mistuning) gives every blade a slightly different
  `f0` and damping. Each blade also sits at a slightly different static
  angular position (stagger error), and the blades untwist slowly with speed.

## Excitation

Stationary flow distortion from the upstream vanes and struts excites the
blades at an integer engine order (EO). The build record does not document
which engine order drives the resonance crossed in this run.

## Run

This was a fast acceleration. The shaft runs up from about 6550 rpm to
7550 rpm in about 6 s, and the speed is not perfectly smooth. The blades'
damping is low, typical of integrally bladed rotors.

## Sensors

- **Once-per-revolution (OPR) sensor:** `opr.csv` holds one time stamp per
  revolution, in seconds.
- **Five optical tip-timing probes** are mounted in the casing. Their angular
  positions are measured in the direction of rotation from probe 1:

| probe | file | angle (deg) |
|---|---|---|
| 1 | `probe_1.csv` | 0.0 |
| 2 | `probe_2.csv` | 14.0 |
| 3 | `probe_3.csv` | 48.5 |
| 4 | `probe_4.csv` | 56.5 |
| 5 | `probe_5.csv` | 153.0 |

  The installation tolerance on each probe angle is ±0.1°.
- **Probe files:** each file lists every trigger time (s) for that probe, in
  time order. The trigger electronics occasionally miss a blade pass, and now
  and then they fire on noise. One probe's fibre connection was intermittent
  during this run.
- **Timing:** jitter is about 0.3 µs rms. All channels share one clock.

## Blade numbering

At the instant of the OPR pulse, the tip of **blade 1** is aligned with probe
1 (to within its static stagger). Blade numbers increase in the order the
blades pass any given probe. Blade 2 reaches probe 1 one blade pitch after
blade 1, and so on up to blade 24.
