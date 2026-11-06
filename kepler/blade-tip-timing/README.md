# blade-tip-timing

The deliverable is a reusable analysis program, `/app/btt_analyze.py`, for a
blade-tip-timing test facility. It takes one run directory, which holds the
raw once-per-rev and probe trigger times plus a `rig.json` of run constants,
and reports four things:

- the engine order that drives the resonance;
- each blade's resonance frequency;
- each blade's modal damping ratio;
- each blade's steady-state resonant tip amplitude, the input to a fatigue
  assessment.

The program is run on the visible example run and on three held-out runs.
The held-out rotors differ in blade count, tip radius, Southwell coefficient,
probe count and layout, engine order, speed range and sweep direction.

## Difficulty

The data are synthetic. They come from the disclosed, seeded generator in
`authoring/provenance/generate.py`. Each blade is an uncoupled
single-degree-of-freedom mode with its own natural frequency and damping
(mistuning). The mode follows a Southwell stiffening law, and engine-order
forcing drives it. The generator integrates the response in time through a
fast speed sweep with small fluctuations. Probe arrival times come from the
deflected tip positions. On top of that, the generator adds:

- timing jitter;
- unknown probe-angle installation errors;
- blade stagger and untwist;
- missed and spurious probe triggers;
- a probe dropout across part of the resonance region;
- a faulty once-per-rev record, with missed and spurious OPR pulses.

The program has to cover the full tip-timing reduction for any rotor in the
format, not just the one run the agent can see:

- repair the OPR record;
- assign triggers to blades so that bad triggers do no harm;
- zero the static positions;
- fit the synchronous response across the probes;
- identify the undocumented engine order.

Three physics points decide the numbers:

1. **Centrifugal stiffening.** A speed bandwidth converts to damping through
   `EO − S/EO`. Ignoring this biases damping by roughly +5% to +25%,
   depending on EO and S.
2. **Transient passage.** The sweeps (about 150–200 rpm/s) are fast compared
   with the blades' decay time, because damping is only 0.12–0.35%. The
   resonance passage is therefore transient: the peak shifts, broadens and
   rings. Steady-state curve fits bias frequency by about 0.3% and damping by
   about +7% to +19% on average. The direction of the distortion flips
   between an acceleration and a deceleration, so the program needs a
   transient model driven by the measured speed history.
3. **Attenuated peaks.** The measured peaks are too low, by roughly 5–16% on
   average and up to about 34% for the lightest-damped blades. The
   steady-state resonant amplitude has to be reconstructed from the
   identified forcing gain and damping.

## Reference solution

`solution/solve.py`. `solution/solve.sh` installs it as
`/app/btt_analyze.py` and runs it on the example. The steps:

1. **Read the rig constants** from `rig.json`.
2. **Repair the OPR record** against the local median period: spurious
   mid-revolution triggers are dropped and missed pulses are filled in.
3. **Shaft angle.** A cubic spline through the repaired OPR gives the shaft
   angle at every trigger.
4. **Blade assignment.** Each trigger goes to the nearest blade slot.
   Triggers off-slot and duplicate triggers are dropped.
5. **Deflection.** Tip deflection is the tip radius times the angular lag
   behind the slot. A constant per blade and probe is removed.
6. **Per-revolution fit.** For each revolution and blade, fit
   `c + a·cos(EO·α_p) + b·sin(EO·α_p)` across the probes. The engine order is
   the integer in 1–30 that leaves the least residual in the vibrating
   revolutions. Response jumps are cleaned with a rolling median.
7. **Starting point.** A steady-state, centrifugally stiffened SDOF fit gives
   initial values.
8. **Transient fit.** For each blade, the slow complex envelope of the SDOF
   response is integrated under the measured speed history and fitted to every
   probe sample.
   - Fitted nonlinearly: `f0` and `ζ`.
   - Solved linearly (variable projection): the complex gain, per-probe
     offsets and a quadratic drift.
   - Samples beyond 5σ are trimmed as outliers.

   The outputs are `f_res = EO·f0/√(EO² − S)` and the steady resonant
   amplitude `|G|/(2ζ)`.

Each run takes about 3–10 s.

## Verification

`tests/test.sh` runs pytest in the separate, offline verifier image. Each run
directory under `tests/heldout/` is copied to a scratch directory: the
example plus the three unseen rotors. `python3 /app/btt_analyze.py <run_dir>
<out.json>` is then executed on each copy as the unprivileged user `nobody`,
with a 10-minute limit per run. Ground truth sits in `tests/truth/`, which is
root-only, so the program under test cannot read it.

Every run must meet all of the following:

- the exact engine order;
- every blade's frequency within 0.25% relative;
- every blade's damping within 25% relative;
- the run's mean signed relative damping error within ±7%;
- every blade's steady-state resonant amplitude within 15% relative.

Non-numeric, non-finite or non-positive values and missing blades fail.

**Where the truth comes from.** Ground truth is the generator's own
parameters, independent of the reference solution:

- the engine order;
- `f_res = EO·f0/√(EO² − S)`;
- `ζ`;
- the forcing deflection divided by `2ζ`.

**Tolerance study.** `authoring/provenance/check.py` and `seed_sweep.sh`
reproduce it.
- **Reference:** on the 4 shipped runs and 12 re-seeded copies, it gets the
  engine order right every time. Frequency is within 0.05%, damping within
  17.7% on the worst blade (run means within ±1.9%), and amplitude within
  9.6%.
- **Steady-state fit with measured peaks:** fails on every run (frequency
  about 0.3%, amplitude 22–34%).
- **Ignoring Southwell:** fails on every run except the high-EO, low-S
  rotor_c, so it never passes all four runs.
- **Truth-reading program:** a program that tries to read the truth files
  finds none and fails.

A constant answer fails the per-blade checks.
