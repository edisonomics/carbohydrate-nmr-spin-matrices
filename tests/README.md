# Sucrose regression tests

Run the fast repository checks from the repository root:

```bash
python3 -m unittest discover -s tests -v
```

For the full physics test, use the project environment. The current lab
environment is `/Users/cece/miniconda3/envs/sucrose`; a reproducible package
specification is [environment.sucrose.yml](../environment.sucrose.yml). To
make a separate copy without changing the working environment:

```bash
conda env create -f environment.sucrose.yml --name sucrose_project
conda run -n sucrose_project python -m unittest discover -s tests -v
```

To use the existing environment directly:

```bash
conda run -n sucrose python -m unittest discover -s tests -v
```

These tests do not require MATLAB or Spinach. They check the molecule
configuration, the downloaded BMRB/GISSMO provenance, matrix dimensions and
numeric equivalence, the Bubb-seed parser, and the independent simulator's
physics self-tests. A failing matrix/provenance test means the data-selection
step needs attention; it does not by itself mean the experimental sample must
be reacquired.

The identity-free mystery-sugar screen is tested in
`test_mystery_sugar_workflow.py`. It verifies that a two-anomer 1-D fingerprint
can rank a xylose-like candidate while a reference-free fructose family stays
in `REVIEW`.

`test_carbohydrate_input_workflow.py` tests the reviewed four-file bootstrap:
one-to-one atom/spin mapping, component and coupling references, same-sample
multifield design, held-out validation, deterministic symmetric matrices, and
the rule that existing provisional output is not overwritten silently.
