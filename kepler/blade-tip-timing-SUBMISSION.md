# blade-tip-timing: submission notes

## Form fields

- **Domain:** Engineering Sciences
- **Field:** Mechanical & Aerospace Engineering
- **GPU:** No (CPU only)
- **Task name:** `blade-tip-timing`
- **Bundle:** `kepler/blade-tip-timing.zip`

## Phase-1 idea check (paste text)

A compressor blisk with 24 mistuned blades is accelerated quickly
(6550→7550 rpm in about 6 s) through a first-bending resonance. The only data
are raw blade-tip-timing records: once-per-rev time stamps and the trigger
times from five casing probes. The solver must find the exciting engine order
(undocumented), plus each blade's resonance frequency and modal damping ratio.

The core work is real BTT practice:

- turn trigger times into tip deflections;
- assign pulses to blades despite missed and spurious triggers and a probe
  that drops out mid-resonance;
- remove static stagger and probe-angle offsets;
- run a per-revolution multi-probe synchronous fit and identify the engine
  order from it;
- fit a damped single-mode model per blade.

The difficulty is physical, not puzzle-like:

1. The mode stiffens centrifugally (Southwell coefficient given). The speed
   bandwidth then maps to damping through `EO − S/EO`, not `EO`. Ignoring
   this biases damping by about +13%.
2. The run-up is fast compared with the blades' decay time (damping
   0.12–0.35%). The resonance passage is therefore transient: the peak
   shifts, broadens and rings afterwards. A steady-state fit or half-power
   method biases frequency by about 0.3% and damping by about +15%. Damping
   has to come from a transient model driven by the measured speed history.

Requirements (stated in the instruction):

- correct engine order;
- every blade's frequency within 0.25%;
- every blade's damping within 25%;
- damping error averaged over the blades within ±7%.

The reference transient fit lands within 0.04% on frequency and within 9% on
damping per blade (mean within ±1.6%). This holds on the shipped data and on
3 freshly generated datasets.

## Self-review against the posted rejection reasons

| Rejection reason | Status |
|---|---|
| Tolerance inside the noise of correct methods | Reference worst blade 7–9% against a 25% limit. Mean ±1.6% against ±7%. Frequency 0.04% against 0.25%. Checked on 4 datasets. |
| Accuracy bars not visible to the solver | All four limits are stated in instruction.md. |
| Unstated conventions | Blade numbering, probe angles and rotation direction, the damping definition (viscous single-mode, fraction of critical), the frequency definition (natural frequency at the crossing speed) and the Southwell law are all stated in the rig sheet and instruction. |
| Wrong shipped data without a warning | None. Probe-angle tolerance (±0.1°), stagger, missed and spurious triggers, and the intermittent probe are all disclosed. |
| Shortcuts that leak the answer | No fitted values ship. The engine order and damping can only come from the probe data. |
| Results don't separate good from bad | Steady-state fit and no-Southwell variants fail on every dataset tried. The transient reference passes on every one. |
| Puzzle or spec transcription | No. It is a measurement-data inference pipeline, and the method is not prescribed. |
| Not a scientific result | Modal properties of an engineered structure. |

## Instruction wording (new rules)

- No grader, verifier, scoring or reward language.
- Headings are `## Output` and `## Requirements`.
- `[agent] timeout_sec = 28800` and `[verifier] timeout_sec = 900`.
- The last line of the instruction states the 8-hour limit.

## Still to confirm before upload

- The `domain` and `field` slugs in task.toml (`engineering-sciences`,
  `mechanical-aerospace-engineering`) must match the form.
- The required-layout line on the form was cut off at `tests/Do…`. If it
  names an extra file, it needs adding.
- The Docker build has not been run. The solution was run with the
  Dockerfile's pinned package versions in a local virtualenv.
