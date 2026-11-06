A small solid rocket motor was static-fired once on a thrust stand. The record is in `/app/data/static_fire.csv` (columns `time_s`, `chamber_pressure_MPa_gauge`, `thrust_N`, with ignition near `time_s` = 0) and the motor, propellant, nozzle and site data are in `/app/data/motor.json`. The grain is a single cylindrical BATES segment whose bore and both end faces burn. The graphite throat eroded during the firing; its diameter was measured before and after the test, but its history during the burn is unknown. Chamber pressure is reported relative to ambient, and the load cell was not re-zeroed before the firing.

Reconstruct the motor's internal ballistics from this one firing. Use a lumped (0-D) chamber with the burning-rate law $ r = r_4 \, (p_c / 4\,\mathrm{MPa})^{n} $, where $ p_c $ is absolute chamber pressure, and a constant delivered characteristic velocity $ c^* $. The nozzle is conical with the divergence half-angle given in `motor.json`, the divergence factor is $ \lambda = (1 + \cos\alpha)/2 $, the exhaust is an ideal gas with the given constant $ \gamma $, and the exit diameter does not change. Thrust is $ F = \eta \, \lambda \, \dot m \, u_e + (p_e - p_a) A_e $, where $ u_e $ and $ p_e $ follow from isentropic expansion through the current area ratio $ A_e / A_t $ and $ \eta $ is an unknown constant jet efficiency. Both the loaded propellant and the igniter gas leave through the nozzle, and no propellant remains after burnout.

Write `/app/output/ballistics.json` containing exactly these keys, with plain JSON numbers in the units named:

- `c_star_m_per_s`: delivered $ c^* $, within 1.5% of the true value.
- `jet_efficiency`: $ \eta $, within 0.01 of the true value.
- `burn_rate_exponent`: $ n $, within 0.02 of the true value.
- `burn_rate_at_4MPa_mm_per_s`: $ r_4 $ in mm/s, within 2% of the true value.
- `total_impulse_N_s`: total impulse of the firing in N·s, within 1% of the true value.
- `throat_diameter_mm`: an object with exactly the keys `"1.0"`, `"1.5"`, `"2.0"`, `"2.5"` and `"3.0"`, giving the throat diameter in mm at those values of `time_s`, each within 0.05 mm of the true value.

All numbers must be finite.

You have 28800 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
