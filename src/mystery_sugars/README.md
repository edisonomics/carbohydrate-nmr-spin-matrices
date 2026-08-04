# Mystery-sugar 1-D identification project

This project is the identity-free front end for an unknown carbohydrate measured at CCRC on multiple 1-D `¹H` fields (for example 600, 900, and 1100 MHz).

It does not assume the molecule name and it does not invent a final matrix. The screen:

1. reads the already prepared DSS-referenced spectra;
2. detects multiplet centers across the complete prepared proton window;
3. compares both the full fingerprint and anomeric region with BMRB references;
4. where a local matrix exists, exactly simulates the multiplet shapes and all
   scalar couplings at every measured field;
5. evaluates source-traceable Bubb/Duus structural, coupling, and evidence rules;
6. ranks topology hypotheses; and
7. writes `CANDIDATE` or `REVIEW` plus the evidence needed next.

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

The default window is 0.50–10.00 ppm from each DSS-aligned full acquisition,
excluding water at 4.65–4.90 ppm. Use
`--no-physics-model` only when you deliberately want the faster shift-only
screen.

The combined score uses 70% complete chemical-shift fingerprint, 20% exact
matrix-derived multiplet shape, and 10% cited literature guidance. A missing optional
matrix or literature channel receives a neutral 0.5 rather than being confused with
negative evidence. Generated BMRB candidates remain review-gated.

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

The expanded mode adds candidates with assigned proton references, including
entries without an anomeric proton. All newly generated candidates remain
`needs_review`. If one ranks first, the report must remain `REVIEW`; only the
hand-reviewed core library can reach `CANDIDATE` from the 1-D screening stage.

## What the multiplet/J score means

For candidates with a local model configuration, the native Python quantum
engine builds the full Hamiltonian from the matrix shifts and nonzero J
couplings, then forward-simulates the complete spectrum at each field. The
reported cosine similarity therefore tests the multiplet shapes and all J
values together. Linewidth and a small global ppm offset are nuisance searches.

This score is screening evidence. It does not assign individual lines or turn
an unresolved spacing into a signed scalar coupling. Candidates without a
local matrix show `n/a` rather than receiving invented coupling evidence.

## Executable literature knowledge

The machine-readable source of truth is
`src/common/carbohydrate_nmr_knowledge.json`. It currently contains rules
verified directly against Bubb (2003), DOI `10.1002/cmr.a.10080`, and Duus,
Gotfredsen & Bock (2000), DOI `10.1021/cr990302n`. Every rule records:

- a stable rule ID;
- the paper DOI and exact section/page locator;
- the scientific statement and its applicability conditions;
- quoted numeric values separately from software tolerances; and
- the permitted use in scoring, warnings, or confirmation gates.

The current executable rules report:

- resolved anomeric reporters and the visible component count;
- methyl reporters near 1.2 ppm and acetyl reporters near 2.0–2.1 ppm;
- cross-field anomeric line spacings;
- general alpha-like 2–4 Hz and beta-like 7–9 Hz J1,2 patterns; and
- the mannose exception, approximately 1.6 Hz for alpha and 0.8 Hz for beta.

Only spacings reproduced in at least two fields enter this guidance channel.
Bubb explicitly warns that crowded nonanomeric carbohydrate spectra are
often not first order, so ordinary line separations must not automatically be
called coupling constants. Duus rules add stereochemical coupling checks,
optional one-bond carbon-proton checks, sample-metadata requirements, and the
need for multidimensional confirmation. Use `--no-literature-guidance` to
disable this channel for a controlled comparison; `--no-bubb-guidance` remains
as a backward-compatible alias.

The JSON report retains every rule result with its observation, score, status,
source, locator, applicability, and explanation. Adding another paper means
adding verified rules to the knowledge file and regression tests—not embedding
uncited prose in the ranker.

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
