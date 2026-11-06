I need creep-force curves for a worn wheel running on our rail section, computed with Kalker's exact rolling-contact theory: the elastic half-space, boundary-element approach that tools like CONTACT use. Simplified models such as Hertz plus FASTSIM aren't accurate enough for this profile pair.

The inputs are in `/app/data/`:

- `/app/data/rail_profile.csv` and `/app/data/wheel_profile.csv`: the measured transverse profiles of the rail head and the wheel tread, in the same frame. The columns are lateral position y and surface height z, both in mm.
- `/app/data/contact_case.json`: materials, normal load, friction coefficient, wheel rolling radius, the frame and sign conventions for creepage and force, and the eight creepage cases A to H, each with longitudinal ξ, lateral η and spin φ.

Treat both bodies as linearly elastic half-spaces made of the same steel, under steady-state rolling, with Coulomb friction at the given coefficient. The contact is planar.

Write `/app/output/creep_forces.json` in this format:

```json
{
  "contact_area_mm2": 0.0,
  "max_pressure_MPa": 0.0,
  "cases": [
    {"id": "A", "Fx_N": 0.0, "Fy_N": 0.0}
  ]
}
```

`contact_area_mm2` and `max_pressure_MPa` describe the normal contact under the given load. `cases` holds one entry for each of the eight creepage cases, with the total tangential force that the rail exerts on the wheel, using the sign convention in `contact_case.json`.

Requirements:

- For every case, the force vector (Fx, Fy) is within 2.5% of μN of the exact-theory solution, where μN is the friction coefficient times the normal load.
- The contact area is within 4% of the exact-theory value.
- The maximum contact pressure is within 5% of the exact-theory value.
- All eight cases are present, and every value is a finite number.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
