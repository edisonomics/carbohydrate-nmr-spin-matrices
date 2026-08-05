# Reviewed carbohydrate input workflow

Use this workflow when the carbohydrate is known and you have proton spectra at
multiple magnetic fields, but no verified spin-system matrix. It turns reviewed
assignments and coupling evidence into a **provisional** starting matrix. It does
not infer atom identities from peak order and it never overwrites a matrix under
`data/`.

## The workflow

```text
known structure + one physical sample
                 |
                 v
  molecule.json                  solution forms and conditions
  proton_assignments.csv         one row for every modeled proton
  couplings.csv                  one row for every modeled J pair
  spectra.csv                    training and held-out spectra
                 |
                 v
       strict input validation
                 |
                 v
 provisional matrices + fit configuration + provenance report
                 |
                 v
 manual mapping review -> multifield fit -> held-out-field validation
```

The manual review before fitting is essential. A syntactically valid file can
still attach a correct chemical shift to the wrong proton.

## 1. Create the four input files

From the repository root:

```bash
python3 src/common/bootstrap_carbohydrate.py init \
  --molecule galactose \
  --forms alpha_pyranose beta_pyranose
```

This creates `data/galactose/input/`. That directory is ignored by Git because
it may refer to local experimental data. Use `--input-dir` to put the bundle
somewhere else.

### `molecule.json`

Record the chemical identity, every solution form included in the model, the
fraction of each form, and the experimental conditions. Fractions must sum to
1.0. Use `PROVISIONAL` when a fraction is only an initial estimate.

Do not combine alpha and beta anomers into one component. Ring forms and other
slowly exchanging species that have distinct proton resonances also need
separate components.

Required condition fields are solvent, temperature in kelvin, pH, and the
chemical-shift reference. These are part of the data, not optional notes:
chemical shifts and anomer populations can change with conditions.

### `proton_assignments.csv`

Each row represents exactly one modeled proton.

| Column | Meaning |
| --- | --- |
| `component` | Exact form ID from `molecule.json`. |
| `spin_id` | Unique internal ID, such as `alpha_H1`. It must be unique across the whole bundle. |
| `atom_label` | Human-readable proton name in the chosen structure, such as `H1` or `H6a`. |
| `attached_atom` | The heavy atom carrying that proton, such as `C1` or `C6`. |
| `bmrb_label` | Original deposited atom label, if a BMRB assignment is used. Preserve it rather than silently renaming it. |
| `shift_ppm` | Initial proton chemical shift under comparable conditions. |
| `assignment_status` | Evidence category defined below. |
| `source` | Entry ID, DOI/table, spectrum/peak ID, or an explicit lab assignment note. |
| `uncertainty_ppm` | Honest uncertainty on the starting value. |
| `fit` | `true` if the multifield optimizer may refine the shift. |

For methylene groups, enter `H6a` and `H6b` as separate spins. Do not decide
which is pro-R/pro-S merely from left-to-right peak order. If stereospecific
assignment is unavailable, keep a documented `a/b` convention, use a status
that reflects the uncertainty, and review both mappings during validation.

### `couplings.csv`

Each row represents one scalar coupling in one component. `spin_i` and
`spin_j` must refer to spin IDs in `proton_assignments.csv`. The generated
matrix is symmetric, so list each unordered pair only once.

`J_hz` must come from an assigned multiplet, a J-resolved experiment, a
deposition that explicitly reports J, or a cited literature prior. A 2-D
cross-peak supports connectivity but does not by itself measure J. Likewise,
the distance between two unrelated lines is not a coupling.

Use a signed value only when the sign convention is supported by the source.
Record the uncertainty, evidence category, source, and whether the optimizer
may refine the value.

### `spectra.csv`

List the one-dimensional proton spectra used in the joint fit.

- Supply at least two distinct `training` fields and one independent
  `validation` field.
- All rows must have the same `sample_id` and `tube_id`. This prevents
  concentration, pH, temperature, or composition differences from being
  mistaken for magnetic-field effects.
- The validation field cannot also be a training field.
- `path` may be relative to the repository root or absolute. It must point to a
  processed Bruker experiment directory containing `acqus` and
  `pdata/<procno>/procs` plus `pdata/<procno>/1r`. Record that processed-data
  number in the `procno` column. Use `import_bruker_dataset.py` first if the
  experiment has not yet been copied into this repository.
- Record the actual pulse program; do not label an edited or decoupled
  experiment as an ordinary 1-D proton spectrum.

## Evidence categories

| Status | Use |
| --- | --- |
| `DEPOSITED` | Directly reported in a curated deposition for this molecule/form. |
| `ASSIGNED` | Assigned from this sample with traceable peak or experiment evidence. |
| `LITERATURE_PRIOR` | Reported for this molecule/form in literature, but not independently measured in this sample. |
| `ANALOG_PRIOR` | Borrowed from a related molecule; always produces a review warning. |
| `UNKNOWN` | Placeholder with no adequate evidence; always produces a review warning. |

Form fractions use `MEASURED`, `ASSIGNED`, `PROVISIONAL`, or `UNKNOWN`.

## 2. Validate before generating anything

```bash
python3 src/common/bootstrap_carbohydrate.py build \
  --input-dir data/galactose/input \
  --validate-only
```

The validator checks names, one-to-one spin IDs, component membership, duplicate
couplings, numeric ranges, fraction totals, source fields, same-sample design,
required Bruker files, the stated field/nucleus/pulse program against `acqus`,
and the training/validation split. Its report is
written to `outputs/galactose/bootstrap/`.

`PASS` means only that the bundle is structurally consistent enough to build a
provisional seed. It does **not** prove the molecular identity, atom mapping,
coupling sign, or assignments.

## 3. Generate provisional artifacts

After resolving all errors and reviewing warnings:

```bash
python3 src/common/bootstrap_carbohydrate.py build \
  --input-dir data/galactose/input
```

The command writes:

```text
outputs/galactose/bootstrap/
|-- matrices/
|   |-- alpha_pyranose_provisional_matrix.txt
|   `-- beta_pyranose_provisional_matrix.txt
|-- galactose_config.provisional.json
|-- fit_parameter_policy.json
|-- input_validation_report.json
|-- input_validation_report.md
`-- seed_manifest.json
```

The diagonal of each matrix contains chemical shifts in ppm. Symmetric
off-diagonal entries contain J couplings in hertz. Connected spin blocks are
generated from the coupling graph using one-based indices expected by the
current simulator configuration.

The command refuses to replace an existing bootstrap output unless `--force`
is supplied. `--force` is for regenerating reviewed provisional output; it
still cannot write into the canonical `data/<molecule>/` matrix location.

## 4. Scientific review and fitting

Before fitting, print or inspect a mapping table containing component, spin ID,
atom label, attached carbon, BMRB label, shift, coupling partners, and source.
Check it against the molecular drawing and any HSQC/COSY assignments. The most
important checks are:

1. every physical proton appears exactly once in the intended component;
2. methylene partners and anomers are not silently exchanged;
3. shifts match the stated solvent, pH, temperature, and reference closely
   enough to be valid starting values;
4. every nonzero matrix entry has traceable J evidence;
5. only training fields are optimized;
6. the untouched validation field improves after refinement, including full
   spectrum, anomeric region, multiplet shapes, and residuals.

Only after those checks should a provisional matrix be promoted into a
molecule-specific fitting workflow. GISSMO submission or replacement of a
verified matrix remains a separate, manual decision.
