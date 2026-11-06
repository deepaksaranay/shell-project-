# blade-tip-timing

A mistuned 24-blade compressor rotor is run up quickly through a first-bending
resonance while five casing probes time every blade-tip pass. From the raw
once-per-rev and trigger times, the task is to identify three things: the
exciting engine order, and each blade's resonance frequency and modal damping
ratio.

## Difficulty

The data are synthetic. They come from the disclosed, seeded generator in
`authoring/provenance/generate.py`. Each blade is an uncoupled
single-degree-of-freedom mode with its own natural frequency and damping
(mistuning). The mode follows a Southwell stiffening law, and engine-order
forcing drives it. The generator integrates the response in time through a
speed history with small fluctuations. Probe arrival times come from the
deflected tip positions. On top of that, the generator adds:

- timing jitter;
- unknown probe-angle installation errors;
- blade stagger and untwist;
- missed and spurious triggers;
- a probe dropout over part of the resonance region.

The agent has to build the full tip-timing data reduction. That covers blade
assignment that survives bad triggers, static zeroing, the multi-probe
synchronous fit, and identifying the undocumented engine order. Two physics
points decide whether the damping comes out right:

1. **Centrifugal stiffening.** The mode stiffens with speed, so a speed
   bandwidth converts to damping through `EO − S/EO`. Ignoring this gives
   about +13% bias.
2. **Transient passage.** The run-up (about 167 rpm/s) is fast compared with
   the blades' decay time, because damping is only 0.12–0.35%. The resonance
   passage is therefore transient: the peak shifts, broadens and rings
   afterwards. Steady-state curve fitting and the half-power method bias the
   frequency by about 0.3% and the damping by about +15%. Damping has to come
   from a transient model driven by the measured speed history.

## Reference solution

`solution/solve.py`, called from `solution/solve.sh`:

1. **Shaft angle.** A cubic spline through the OPR times gives the shaft
   angle at every trigger.
2. **Blade assignment.** Each trigger goes to the nearest blade slot.
   Triggers far from any slot, and duplicate triggers, are dropped.
3. **Deflection.** Tip deflection is the tip radius times the angular lag
   behind the slot. A constant per blade and probe is removed.
4. **Per-revolution fit.** For each revolution and blade, fit
   `c + a·cos(EO·α_p) + b·sin(EO·α_p)` across the probes. The engine order is
   the integer that leaves the least residual in the vibrating revolutions.
5. **Starting point.** A steady-state, centrifugally stiffened SDOF fit gives
   initial values for the final step.
6. **Transient fit.** For each blade, the slow complex envelope of the SDOF
   response is integrated under the measured speed history and fitted to every
   probe sample.
   - Fitted nonlinearly: `f0` and `ζ`.
   - Solved linearly (variable projection): the complex gain, per-probe
     offsets and a quadratic drift.
   - Samples beyond 5σ are trimmed as outliers.

   The output is `f_res = EO·f0/√(EO² − S)`.

The full run takes about 5 s.

## Verification

Ground truth is the generator's own parameters, independent of the reference
solution. It is stored in `tests/expected.json`, which is baked into the
verifier image. The tests read `/app/output/blade_modal_results.json` and
require all of the following:

- the exact engine order;
- every blade's frequency within 0.25% relative;
- every blade's damping within 25% relative;
- the signed relative damping error averaged over the 24 blades within ±7%.

Non-numeric, non-finite or non-positive values fail.

`authoring/provenance/check.py` and `seed_sweep.sh` reproduce the tolerance
study below. Each row covers the shipped data plus 3 freshly generated seeds.

| Method | Frequency error | Damping, worst blade | Damping, mean bias | Result |
|---|---|---|---|---|
| Transient reference | ≤ 0.04% | ≤ 9.4% | within ±1.6% | Passes on every dataset |
| Steady-state SDOF fit | ≈ 0.3% | — | +8.5% to +15% | Fails on every dataset |
| Southwell ignored | — | — | +12% to +15% | Fails on every dataset |

A constant answer fails the per-blade checks.
