# surge-vessel-sizing

The task is a hydraulic-transient design problem. A pump in a 5.3 km rising
main trips. The agent has to size an air vessel at the pump station: total
volume, initial gas volume and inflow-orifice diameter. The vessel must stay
within a volume budget, and after the trip it must keep the whole pipeline
within its pressure limits and stop the vessel from draining.

## Difficulty

The system is a realistic but synthetic design case. Its data are in
`environment/data/system.json`. No ready-made simulator or answer ships.
The agent has to build a transient model good enough to design with:

- method of characteristics in three pipes with different diameters and
  wave speeds, plus series junctions;
- quasi-steady friction;
- the steady operating point from the pump curve and the system curve;
- pump run-down from the homologous head and torque curves and the rotor
  inertia;
- an ideal check valve, including its closure logic once the pump can no
  longer beat the line pressure;
- an air-vessel boundary with polytropic gas, the changing water level and
  differential connection losses;
- a check of the pressure envelope along the whole route profile.

Then it has to search the three-dimensional design space against three
interacting constraints:

- **Minimum pressure.** It is set by the mass oscillation over the knee at
  chainage 1500 m, about 30 s after the trip.
- **Gas fraction.** Too much initial gas drains the vessel.
- **Return upsurge.** An orifice larger than about 0.19 m lets the upsurge
  exceed 105 m.

The feasible minimum under the verifier model is about 16.9 m³. The budget
is 18.0 m³, and at that volume only initial gas fractions of about 0.42–0.46
work.

Shortcuts and modelling mistakes lead to designs that fail:

| approach | design it leads to | result |
|---|---|---|
| Check pressure only at the pump | about 12 m³ | fails at the knee |
| Isothermal gas (n = 1) | about 15 m³ | fails |
| Instant pump stop (no inertia) | no feasible design under 30 m³ | — |
| 50% initial gas, the usual rule of thumb | — | the vessel drains past 90% |
| Large orifice | — | upsurge above 105 m |

## Reference solution

`solution/solve.py` is an independent transient model, written separately
from the checking code. It differs from the checking model in three ways:

| aspect | reference model | checking model |
|---|---|---|
| discretisation | space-line interpolation at Courant number ≤ 1 | wave-speed adjustment |
| node solver | scipy `fsolve` at the pump/check-valve/vessel node | Newton iteration |
| rotor speed | Heun integration | trapezoidal integration |

On five test designs the two models agree to within 0.01 m. The reference
fixes the volume at 17.9 m³, scans the initial gas fraction and two orifice
sizes, and writes the design with the largest slack against every limit:

| design choice | value |
|---|---|
| vessel total volume | 17.9 m³ |
| initial gas volume | 7.966 m³ |
| inflow orifice diameter | 0.12 m |

Its predicted performance:

| quantity | value |
|---|---|
| minimum absolute pressure head | 3.90 m |
| maximum gauge pressure head | 102.4 m |
| maximum gas fraction | 0.870 |

The search takes about 1.5 minutes.

## Verification

`tests/test.sh` runs pytest offline in the verifier image. The tests read
`/app/output/surge_vessel_design.json` and check the budget, the ranges and
the number types. They then simulate the design with `tests/hydraulics.py`,
using MOC at dt = 0.005 s with wave-speed adjustments under 0.5%, over the
full 180 s. This model is not the one in the reference solution.

**What the design must meet:**

- minimum absolute pressure head at least 3.0 m at every computational
  section of all pipes;
- maximum gauge pressure head at most 105.0 m;
- gas fraction at most 0.9.

**Where the truth comes from.** The ground truth is this physical
simulation of the submitted design. The model is converged: the metrics
change by less than 0.05 m between dt = 0.02, 0.01 and 0.005 s, which is
small against the roughly 1 m slack available at the budget.

**Validation.** `authoring/evidence/validation.md` records the oracle (1),
the nop (0) and seven lazy or mistaken designs. All of those designs fail
except one inside the narrow feasible window.
