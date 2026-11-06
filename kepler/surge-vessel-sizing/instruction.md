We need surge protection for a pumped rising main. At present a power failure at the pump station would trip the pump, the line would drain towards the high ground, and the pressure would collapse over the knee and the summit of the route. When the water column came back, the upsurge would hit the pump station. The plan is to add an air vessel on the discharge header, just downstream of the pump's check valve, and size it as small as the budget allows.

Everything known about the system is in `/app/data/system.json`. It holds the fluid properties, the reservoir levels, the pump (homologous head and torque curves, rated point, rotor inertia), the check valve, the three pipes (lengths, diameters, wave speeds, friction factors), the route profile, the air-vessel geometry rules and connection losses, the trip event, the design limits, the volume budget, and the modelling basis for the transient analysis.

Three design choices are open: the vessel's total internal volume, the gas volume in it during normal steady operation before the trip, and the diameter of the inflow orifice that throttles water flowing back into the vessel.

Put the design in `/app/output/surge_vessel_design.json`:

```json
{
  "vessel_total_volume_m3": 0.0,
  "initial_gas_volume_m3": 0.0,
  "inflow_orifice_diameter_m": 0.0
}
```

Requirements:

- The total vessel volume is at most 18.0 m³.
- The inflow orifice diameter is between 0.05 m and 0.35 m.
- The initial gas volume is positive and at most 90% of the total vessel volume.
- Over the 180 s after the pump trip, everywhere along all three pipes, the absolute pressure head never falls below 3.0 m and the gauge pressure head never exceeds 105.0 m.
- The gas never occupies more than 90% of the vessel volume.

The transient has to meet these limits when it is computed as described in the modelling basis in `system.json`. That is a converged one-dimensional elastic water-hammer solution: quasi-steady friction, the pump running down on its own inertia and hydraulic torque, the ideal check valve, and the polytropic gas with the vessel's water level and connection losses.

Python 3 with numpy, scipy and pandas is available.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
