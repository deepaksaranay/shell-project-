# blade-tip-timing: submission notes

## Form fields

- **Domain:** Engineering Sciences
- **Field:** Mechanical & Aerospace Engineering
- **GPU:** No (CPU only)
- **Task name:** `blade-tip-timing`
- **Bundle:** `kepler/blade-tip-timing.zip`. It contains task.toml,
  instruction.md, README.md, environment/, solution/, tests/ and authoring/.
- **Deliverable:** the program `/app/btt_analyze.py`, run on held-out
  datasets. This changed from an answer file to address the similarity
  rejection.

## Phase-1 idea check (paste text)

The task is to build a reusable blade-tip-timing analysis program for a
blisk test facility. Given one run's raw once-per-rev and casing-probe
trigger times, plus a rig.json of rotor constants, the program must find
four things:

- the undocumented engine order;
- each mistuned blade's resonance frequency;
- each blade's modal damping;
- each blade's steady-state resonant tip amplitude, the input to a fatigue
  assessment.

The program is run on the visible example and on three unseen rotors. Those
rotors differ in blade count, radius, Southwell coefficient, probe layout,
engine order, speed range and sweep direction (acceleration or deceleration).

It needs a complete, robust tip-timing reduction:

- repair the OPR record, which has missed and spurious pulses;
- assign triggers to blades despite missed and spurious triggers and a probe
  dropout;
- zero the static positions;
- fit the synchronous response across the probes;
- identify the engine order.

It also needs real physics judgement:

- Centrifugal stiffening maps speed bandwidth to damping through `EO − S/EO`.
- The sweeps are fast relative to the blades' decay time, so the resonance
  passage is transient. Its distortion flips between run-up and run-down,
  and steady-state fits fail.
- Measured peaks are attenuated by up to about 34%, so the fatigue amplitude
  has to be reconstructed from the identified forcing gain and damping.

The data come from a disclosed, seeded generator, and the ground truth is the
generator's parameters.

## Self-review against the posted rejection reasons

| Rejection reason | Status |
|---|---|
| Tolerance inside the noise of correct methods | Over 16 runs (4 configurations × 4 seeds), the reference's worst blade is 17.7% on damping against a 25% limit and 9.6% on amplitude against 15%. Run means stay within ±1.9% against ±7%. Frequency is within 0.05% against 0.25%. |
| Accuracy bars not visible to the solver | All four limits are stated in instruction.md. |
| Unstated conventions | Blade numbering, probe angles and rotation direction, the damping definition (viscous single-mode, fraction of critical), the frequency definition (natural frequency at the crossing speed) and the Southwell law are all stated in the rig sheet and instruction. |
| Wrong shipped data without a warning | None. Probe-angle tolerance (±0.1°), stagger, missed and spurious triggers, the intermittent probe and the faulty OPR pickup are all disclosed. |
| Shortcuts that leak the answer | No fitted values ship. Held-out runs are unseen, and the truth directory is unreadable to the program under test. |
| Results don't separate good from bad | Steady-state fit and no-Southwell variants fail on every dataset tried. The transient reference passes on every one. |
| Puzzle or spec transcription | No. It is a measurement-data inference pipeline, and the method is not prescribed. |
| Not a scientific result | Modal properties of an engineered structure. |

## Instruction wording (new rules)

- No grader, verifier, scoring or reward language. Requirements are
  introduced with a plain "Requirements:" lead-in.
- `[agent] timeout_sec = 28800` and `[verifier] timeout_sec = 600`.
- The instruction has no headings and ends with the exact required sentence
  ("You have 28800 seconds...").

## Still to confirm before upload

1. **Personal details in task.toml.** Three fields are still `TODO`:
   `author_organization`, `author_profile` and `relevant_experience`.
   Placeholders fail the author-facts check.
2. **instruction.md wording.** It must be in your own words: the platform
   screens instructions for AI-generated text and rejects flagged ones. Keep
   every fact, path, number and tolerance. Use no Markdown headings. Keep the
   exact closing sentence.
3. **Local harbor run.** If you can, run harbor locally (oracle, nop and
   `harbor check`). Docker was not available where this bundle was built.
