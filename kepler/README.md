# Kepler task: bearing-impedance-identification

- `bearing-impedance-identification/`: the task bundle (task.toml, instruction.md,
  environment/, solution/, tests/).
- `bearing-impedance-identification.zip`: the same bundle zipped at its root, ready to upload.
- `generator/`: NOT part of the bundle.
  - `generate.py` regenerates the data and `tests/expected.json`.
  - `ablate.py` and `alt_method.py` check that the tolerance separates
    correct methods from flawed ones.
