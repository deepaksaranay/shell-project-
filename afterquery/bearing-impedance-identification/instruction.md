# Identify the force coefficients of a fluid-film journal bearing

A floating-stator bearing test rig was used to measure the dynamic force coefficients of a
fluid-film journal bearing. The rotor spins on rigid supports; the bearing stator floats and is
excited by two orthogonal hydraulic shakers driving a multi-tone (multi-sine) force. Two
proximity probes measure the stator-minus-rotor relative displacement and two accelerometers
measure the absolute stator acceleration.

All data is in `/app/data/`:

- `rig.json` - sample rate, record length, rotor speed, stator mass, and channel descriptions/units.
- `bearing_shake_x.csv`, `bearing_shake_y.csv` - the bearing test, with the shaker driven along
  x and along y respectively (rotor spinning, bearing lubricated).
- `tare_shake_x.csv`, `tare_shake_y.csv` - a tare (baseline) test of the same stator, with the
  rotor removed and no lubricant, so only the stator support structure (pitch stabilizers)
  reacts the motion. Assume rotor motion is negligible, so in the tare test the probe signal
  equals the stator motion.

Columns of every CSV: `time_s, force_x_N, force_y_N, accel_x_g, accel_y_g, disp_x_um, disp_y_um`.
Read `rig.json` carefully for units and for which stator mass applies to which test
(the bearing insert is removed for the tare test, so the stator is lighter).

## Physics

For the stator, Newton's second law in the frequency domain is

    F(w) - m_s * A(w) = [ H_bearing(w) + H_support(w) ] * D(w)

where `F` is the shaker force vector, `A` the stator acceleration vector, `D` the relative
displacement vector, `m_s` the stator mass, and `H` are 2x2 complex impedance (dynamic
stiffness) matrices. The support structure acts in parallel with the bearing, and its impedance
is what the tare test measures. The bearing impedance is modeled with frequency-independent
stiffness, damping and added-mass matrices:

    H_bearing(w) = K + i*w*C - w^2 * M

Notes:

- Each shaker test alone does not determine a full 2x2 impedance; the bearing is cross-coupled
  and each test also produces some force/motion in the orthogonal direction.
- The excitation tone frequencies are not given; determine them from the data. The rotor is
  spinning, so synchronous (running-speed and harmonic) rotor vibration and static eccentricity
  appear in the displacement signals; they are not excitation and must not be used.
- The excitation tones are exactly periodic within each record.

## Task

Identify `K` (N/m), `C` (N*s/m) and `M` (kg) of the bearing alone (support removed).

Write `/app/output/bearing_coefficients.json` (create the directory if needed) with exactly
these keys, as plain JSON numbers in SI units:

```json
{
  "frequencies_hz": [ ... sorted list of the excitation frequencies you used, in Hz ... ],
  "Kxx": 0, "Kxy": 0, "Kyx": 0, "Kyy": 0,
  "Cxx": 0, "Cxy": 0, "Cyx": 0, "Cyy": 0,
  "Mxx": 0, "Mxy": 0, "Myx": 0, "Myy": 0
}
```

Index convention: the first letter is the force direction and the second the displacement
direction, i.e. `Kxy` is the x-direction force per unit y-direction displacement
(row x, column y of `K`).

Your result is graded against the true coefficients: each direct coefficient must be within a
few percent, and cross-coupled coefficients within a few percent of the corresponding direct
coefficient scale.
