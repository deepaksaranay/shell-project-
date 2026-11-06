# Floating-housing bearing test rig: data sheet

## Layout

The test bearing is a plain fluid-film journal bearing pressed into a
floating housing. A rotor runs through the bore. The rotor is carried by
its own rolling-element support bearings in two pedestals. The housing
itself does not touch the rotor except through the oil film.

- The housing is held in place by soft pitch-stabilizer flexures that run
  from the housing to the rig base. The flexure stiffness and damping are not
  documented. The dry baseline test below is how this rig characterises them.
- The rotor and its pedestal supports are not rigid. At the excitation
  frequencies the shaft vibrates measurably in response to the bearing force.
- A static load is applied to the housing in −y by a pneumatic loader through
  a soft spring. The loader adds no measurable dynamic stiffness.
- The floating housing mass is **18.5 kg**. This figure includes the
  bearing, the accelerometers, the probe mounts and the housing-side halves of
  both load cells.

## Coordinates and sign convention

- x is horizontal and y is vertical (up). Angles are measured from +x towards +y.
- `d` is the displacement of the housing **relative to the shaft**,
  `d = [d_x, d_y]`.
- The oil film force acting on the housing is `f = −(K d + C ḋ + M d̈)`.
  `K`, `C` and `M` are 2×2 matrices with `K[i][j]` the force in direction i per unit `d_j`.

## Excitation

Two electrodynamic shakers act on the housing along +x and +y through
stingers. In each test one shaker is driven and the other is left undriven,
but its stinger stays attached. The drive is a multisine with tones at
20, 30, …, 320 Hz (31 tones). Every tone completes a whole number of cycles
in the record.

## Channels (all recordings)

Sampling rate 4096 Hz, 6.0 s per record. The CSV columns are:

| column | sensor | conversion |
|---|---|---|
| `time_s` | time | s |
| `Fx_V`, `Fy_V` | piezo load cells between each stinger and the housing. Each reads the force applied **to the housing** along +x / +y. | 1.124 mV/N |
| `ax_V`, `ay_V` | piezo accelerometers on the housing, +x / +y. These read absolute acceleration. | 100 mV/g, g = 9.80665 m/s² |
| `prox45_V`, `prox135_V` | eddy-current probes mounted on the housing at 45° and 135°. Each probe tip sits on that side of the shaft and faces the shaft centre. | output = −7.874 V/mm × gap (plus a constant). The gap grows when the shaft moves away from the probe tip. |
| `keyphasor_V` | once-per-revolution notch probe on the shaft | the pulse drops to about −9 V once per revolution |

The load cells and accelerometers are AC-coupled, so their DC content is
meaningless. The proximity probes see the shaft surface, so they also pick
up slow-roll runout (mechanical and electrical) that turns with the shaft.

## Tests

- `baseline_dry_x.csv`, `baseline_dry_y.csv`: no oil supply, rotor
  stationary, housing lifted clear of the journal on the flexures. No film
  force acts. Driven shaker x and y respectively.
- `C01_x.csv` … `C08_y.csv`: lubricated, rotating, loaded. See
  `test_log.csv` for each condition's nominal speed, static load and files.
  The logged speeds are the operator's set-points. The actual shaft speed can
  differ slightly, and the keyphasor gives the true speed.
