# Mystery-sugar 1-D identification project

This project is the identity-free front end for an unknown carbohydrate measured at CCRC on multiple 1-D `¹H` fields (for example 600, 900, and 1100 MHz).

It does not assume the molecule name and it does not invent a final matrix. The screen:

1. reads the already prepared DSS-referenced spectra;
2. detects reproducible anomeric multiplet clusters;
3. compares them with literature/BMRB candidate fingerprints;
4. ranks topology hypotheses; and
5. writes `CANDIDATE` or `REVIEW` plus the evidence needed next.

## Run it

From the repository root:

```bash
conda run -n sucrose_project python src/mystery_sugars/identify_from_1d.py \
  --molecule mystery_sugar
```

Prepare the spectra first if needed:

```bash
conda run -n sucrose_project python src/common/prepare_carbohydrate_spectra.py \
  --molecule mystery_sugar
```

Outputs are written to:

```text
outputs/mystery_sugar/identification/identification_1d_report.json
outputs/mystery_sugar/identification/identification_1d_report.md
```

## Interpretation

The candidate library includes xylose, glucose, mannose, sucrose, and a reference-free fructose family. Fructose is deliberately not given a fake numeric fingerprint: the literature shows that it can contain multiple pyranose, furanose, and keto tautomers whose populations depend on conditions.

The screen is a hypothesis ranker. It must not promote an unknown to a publishable identity without a verified BMRB/GISSMO matrix, local COSY/TOCSY/HSQC/HMBC evidence, or a matched authentic standard.

## Expand screening with explore_BMRB

Sync the structure-deduplicated carbohydrate catalog from the sibling
`explore_BMRB` repository (or its live GitHub Pages copy):

```bash
python3 src/mystery_sugars/sync_explore_bmrb_catalog.py
```

Then include its candidates during ranking:

```bash
conda run -n sucrose_project python src/mystery_sugars/identify_from_1d.py \
  --molecule mystery_sugar \
  --include-bmrb-catalog
```

The expanded mode adds candidates with assigned anomeric-region references,
but all newly generated candidates remain `needs_review`. If one ranks first,
the report must remain `REVIEW`; only the hand-reviewed core library can reach
`CANDIDATE` from the 1-D screening stage.

## Georgia-style Mystery Sugar 1 figures

For the saved Mystery Sugar 1/D-xylose candidate fit, generate high-contrast
full-spectrum and multiplet-zoom figures with the Georgia palette:

```bash
conda run -n sucrose_project \
  python src/common/make_mystery_sugar1_georgia_plots.py
```

The figures are written to `outputs/mystery_sugar/georgia_plots/`. They show
the DSS-aligned experiment, the refined alpha/beta Spinach candidate, and the
equal-response candidate baseline. Because xylose has no deposited GISSMO
matrix, the third trace is not labelled GISSMO.
