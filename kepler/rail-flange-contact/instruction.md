I need the contact forces on a flanging wheel in a tight curve: the wheel runs with its flange against the gauge corner of the rail, so it touches the rail in two places, on the tread and on the flange. For each load case I know the vertical wheel load Q and the lateral force Y with which the wheelset presses this wheel outward against the rail, plus the yaw angle, the speed and the wheel's angular velocity. Please work out where the wheel ends up and what each contact patch carries, using Kalker's exact rolling-contact theory (half-space boundary elements) on each patch. Hertz/FASTSIM isn't good enough here.

The inputs are in `/app/data/`:

- `/app/data/rail_profile.csv`: the rail head profile, lateral position y and height z in mm.
- `/app/data/wheel_profile.csv`: the worn wheel profile with its flange, in wheel coordinates, in mm.
- `/app/data/case.json`: the material, the friction coefficient, the tape-circle radius, the frames and sign conventions, how the contact points and patches are defined, how each patch's creepages follow from the wheelset kinematics, the equilibrium conditions, a starting guess for the wheel position, and the six load cases A to F.

The wheel's lateral shift and vertical position are unknown. They follow from equilibrium of the forces the rail exerts on the wheel, with the creep forces included, as stated in `case.json`. Treat both bodies as the same linear elastic steel, in steady rolling with Coulomb friction. Each patch is planar in its own tangent plane, and the patches are elastically independent.

Write `/app/output/flange_contact.json` in this format:

```json
{
  "cases": [
    {"id": "A", "lateral_shift_mm": 0.0,
     "tread": {"contact_angle_deg": 0.0, "normal_force_N": 0.0, "Fx_N": 0.0, "Fs_N": 0.0},
     "flange": {"contact_angle_deg": 0.0, "normal_force_N": 0.0, "Fx_N": 0.0, "Fs_N": 0.0}}
  ]
}
```

`tread` is the patch with the small contact angle and `flange` is the steep one. Use `null` for a patch that is not in contact. `Fx_N` and `Fs_N` are the tangential force components that the rail exerts on the wheel in the patch's local x and s axes, as defined in `case.json`.

Requirements, for every load case:

- The lateral shift is within 0.01 mm of the converged solution.
- Each patch's contact angle is within 0.3°.
- Each patch's normal force is within 1% of Q.
- Each patch's tangential force vector (Fx, Fs) is within 1% of Q.
- The file is a regular JSON file with exactly the keys shown above: top-level `cases`; `id`, `lateral_shift_mm`, `tread` and `flange` in each case; and `contact_angle_deg`, `normal_force_N`, `Fx_N` and `Fs_N` in each patch. Every load case appears once, and every value is a finite number or `null`.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
