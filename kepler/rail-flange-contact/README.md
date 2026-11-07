# rail-flange-contact

The task is a flanging railway wheel in a tight curve. The wheel touches the
rail on the tread and on the flange flank against the rail's gauge corner.
Each of six load cases gives:

- the vertical wheel load Q and the lateral force Y pressing the flange
  against the rail;
- the yaw angle, the speed and the wheel's angular velocity (two cases
  include braking).

From these the task has to find:

- the wheel's equilibrium position (lateral shift and vertical position);
- each contact patch's contact angle, normal force and local creep forces.

Kalker's exact rolling-contact theory is used on each patch.

## Difficulty

**Who does this work.** In practice this analysis is done by railway
vehicle-dynamics engineers and wheel–rail interface specialists at
rolling-stock manufacturers, infrastructure managers and research
institutes. It is used to assess flange contact in tight curves for:

- derailment safety (flange-climb ratios);
- flange wear and gauge-corner damage;
- the track forces fed into multibody simulations.

**Why it is hard.** This is a coupled multi-point contact problem, not one
patch with given creepages. The solver has to:

1. Find the tread and flange contact points from the measured profiles,
   using fixed wheel regions with the separation minimum in each.
2. Set up each patch in its own tangent plane. At the flange that means
   rotating both profiles by a 68° contact angle.
3. Solve displacement-controlled exact normal problems.
4. Derive each patch's creepages from the wheelset kinematics. The flange
   rolls on a radius 9 mm larger than the tread, so it slips at about −2%
   with a spin near 2 /m, while the tread is in partial slip.
5. Solve steady rolling on both patches.
6. Find the unknown wheel position from the equilibrium ΣF_z = Q,
   ΣF_y = Y, where the creep forces feed back into how the load is shared
   between the patches.

Rigid two-point contact is knife-edge sensitive: the flange gap closes by
about 2.5 mm per mm of lateral shift. Only the elastic equilibrium gives a
well-defined answer.

The shortcuts fail (scored through the verifier; details in
`authoring/evidence/validation.md`):

| shortcut | what goes wrong | verifier |
|---|---|---|
| creep forces left out of the equilibrium | flange normal force off by up to 19 kN (24% of Q) | 0/6 |
| one rolling radius for both patches | flange longitudinal force changes sign | 0/6 |
| horizontal patches (no rotation into the tangent plane) | shift off by 0.03 mm, creep forces off by up to 1.8 kN | 1/6 |
| exact theory on a 1 mm grid | flange creep force misses | 4/6 |

The profiles and loads are synthetic but realistic (worn tread, R13 gauge
corner, 68° flank, Q = 60–80 kN, Y/Q = 0.25–0.6), so no published result
applies.

## Reference solution

`solution/solve.py` takes about 70 s for all six cases:

1. **Contact search.** In each wheel region (tread y_w < 12 mm, flange
   y_w ≥ 12 mm) it finds the minimum of the vertical separation on a
   0.01 mm resampling, and takes the contact angle from the rail slope
   there.
2. **Patch gap.** It rotates both profiles into the patch's tangent plane
   and adds the longitudinal curvature cos γ / r_i.
3. **Normal problem.** It uses 0.4 mm half-space boundary elements, with
   closed-form Boussinesq/Cerruti rectangle integrals looked up from offset
   tables, and a displacement-controlled active set.
4. **Tangential problem.** It computes the creepages from the kinematic
   formulas and finds steady rolling as the fixed point of Kalker's
   transient scheme, with FISTA steps warm-started between evaluations.
5. **Equilibrium.** A Newton stage on a normal-force-only model gets close,
   then Broyden's method on the full model (creep forces included) finishes
   the solve.

The patch solver comes from the earlier single-patch version of this
problem. There it was validated against Kalker's tabulated creep
coefficients for a circular contact (`authoring/provenance/validate_hertz.py`).
With ν = 0.25 it gives C11 = 4.17 against the table's 4.12, and
C22 = 3.73 against 3.67.

Through the verifier it passes 6/6 and uses at most 29% of any limit:

| check | worst error |
|---|---|
| normal force | 165 N |
| creep force | 139 N |
| contact angle | 0.088° |
| lateral shift | 0.0004 mm |

## Verification

`tests/test_outputs.py` reads `/app/output/flange_contact.json`. It must be a
regular file (no symlinks) with exactly the documented keys and every case
once. A patch that is not in contact is `null`. For each load case:

| quantity | limit |
|---|---|
| lateral shift | 0.01 mm |
| contact angle (each patch) | 0.3° |
| normal force (each patch) | 1% of Q |
| tangential force vector (Fx, Fs) (each patch) | 1% of Q |

**Where the answer key comes from.** `tests/answer_key.json` comes from
`authoring/provenance/twopoint.py` on 0.3 mm elements. The generator uses:

- parabolic-refined contact points;
- dense influence matrices;
- an active-set normal solve;
- FISTA transient-scheme steps;
- an equilibrium iteration with frozen creep-to-normal force ratios.

Its 0.4 mm run passes using at most 23% of any limit. The shipped solution
doesn't produce the key.

**How the limits were set.** The limits come from how closely structurally
different correct implementations agree. Each one was scored through the
verifier:

| implementation | how it differs | largest share of a limit | result |
|---|---|---|---|
| `authoring/provenance/independent_submission.py` | separately written; bounded-minimisation contact points, non-square 0.25 × 0.35 mm offset elements, FFT convolution, Polonsky–Keer CG normal solver, ADMM transient steps, alternating equilibrium/creep iteration | 15% | 6/6 |
| `solution/solve.py` (shipped) | 0.4 mm, active-set normal, FISTA, Broyden equilibrium | 29% | 6/6 |
| generator at 0.4 mm | coarser grid | 23% | 6/6 |

The independent submission's worst errors are 85 N in normal force, 105 N
in creep force, 0.045° in contact angle and 0.0002 mm in shift.

Correct methods sit well inside the limits, and the shortcuts above fail.
A missing file, a symlink, a copy of the answer key, or extra or renamed
keys all score 0.

Outputs are in `authoring/evidence/`.
