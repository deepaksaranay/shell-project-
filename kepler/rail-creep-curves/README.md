# rail-creep-curves

The task is to compute wheel–rail creep forces for a measured worn-wheel /
rail-head profile pair with Kalker's exact rolling-contact theory, the
half-space boundary-element method. It runs over eight creepage cases that
combine longitudinal, lateral and spin creepage, and also reports the normal
contact area and the peak pressure.

## Difficulty

**Who does this work.** In practice this calculation is done by railway
vehicle-dynamics engineers and wheel–rail interface specialists at
rolling-stock manufacturers, infrastructure managers and railway research
institutes. When a fleet's wheels wear hollow, they compute exact
creep-force curves for the measured worn profile on the rail section in
service. Those curves feed multibody vehicle simulations that assess:

- curving behaviour;
- hunting stability;
- wear and rolling-contact-fatigue risk;
- when to re-profile the wheels.

Using Kalker's exact theory instead of Hertz/FASTSIM is the expert step for
conformal, non-elliptic worn contacts like this one.

The profiles and operating data are synthetic but realistic. The rail crown
is made of R300 and R80 arcs, and the hollow-worn wheel tread has a concave
R360 band plus a shallow wear groove. They ship in
`environment/data/{rail_profile,wheel_profile}.csv` and
`environment/data/contact_case.json`.

The worn tread is nearly conformal with the rail crown. That gives a wide,
flat-topped, asymmetric patch of about 206 mm² and 22 × 11 mm, with a peak
pressure of about 816 MPa. Because the patch is non-elliptic, the usual
engineering shortcuts fail:

| shortcut | contact area | peak pressure | worst creep-force error (share of μN) |
|---|---|---|---|
| Hertz ellipse from the local curvature | 11.6% high | 10% low | 3.6% |
| Kalker linear theory (no slip region) | — | — | up to 45% |
| Exact method on a 1 mm grid | 4.8% high | 4.4% high | 3.7% |

The agent has to implement and converge the exact theory itself. The steps:

- build the gap from the two profiles;
- solve the normal problem with half-space influence coefficients;
- solve the steady-rolling tangential problem, covering stick/slip, the
  coupled lateral and spin terms and the transport derivative.

The naive iterative solver diverges for this non-symmetric problem: block
Gauss–Seidel has a spectral radius of about 20–100. So a robust formulation
is needed, and so is a grid fine enough for about 1% accuracy.

## Reference solution

`solution/solve.py` takes these steps:

1. Builds the undeformed separation.
2. Solves the normal problem on 0.4 mm rectangular elements, using the
   closed-form Boussinesq rectangle integrals and an active-set method with
   the load as a side condition.
3. Gets steady rolling as the fixed point of Kalker's transient rolling
   scheme, with a rolling step of one element.

   Each step is the convex problem: minimise ½ qᵀUq + qᵀ(w·Δx − u_prev) such
   that |q| ≤ μp. U is the symmetric Cerruti influence matrix, and the step
   is solved by accelerated projected gradient (FISTA). Steps repeat until
   the tractions stop changing.

Run through the verifier, it scores as follows:

| check | result |
|---|---|
| worst creep force | 1.81% of μN (pure-spin case F) |
| contact area | +1.70% |
| peak pressure | +1.58% |

The run takes about 11 minutes.

The same scheme is validated in `authoring/provenance/validate_hertz.py`
against Kalker's tabulated creep coefficients for a circular contact:

| ν | quantity | computed | Kalker table |
|---|---|---|---|
| 0 | C11 = C22 | 3.45 | 3.40 |
| 0.25 | C11 | 4.17 | 4.12 |
| 0.25 | C22 | 3.73 | 3.67 |

It also saturates correctly at μN.

## Verification

`tests/test_outputs.py` reads `/app/output/creep_forces.json`. Every value
must be numeric and finite, and every case must be present exactly once.

| quantity | limit |
|---|---|
| each case's creep-force vector | within 2.5% of μN (962 N) of the reference |
| contact area | within 4% |
| peak pressure | within 5% |

**Where the reference comes from.** `tests/reference_solution.json` holds
the exact solution computed with `authoring/provenance/problem.py` at
element sizes of 0.5 mm and 0.35 mm, then Richardson-extrapolated (first
order in element size). The reference solution doesn't produce these values.

**Independent implementation.** `authoring/provenance/independent_check.py`
uses a different tangential algorithm: it solves the steady-state equations
s = w + Kq directly, with an active-set method. It runs on a 0.45 mm grid
shifted by a quarter element. Through the verifier it passes all 10 checks:

| check | result |
|---|---|
| worst creep force | 1.75% of μN |
| contact area | +1.51% |
| peak pressure | +0.70% |

Its output is in `authoring/evidence/independent_output.json`.

**How the limits were set.** The worst force errors against the
extrapolated reference are 1.6% of μN at 0.5 mm and 1.1% at 0.35 mm on the
generator's own grids. The shipped solution is at 1.81% and the independent
implementation at 1.75%. The Hertz-ellipse shortcut fails the area and
pressure limits and several force cases. The 1 mm grid and linear theory
fail the force limits.
