# Sucrose spectrum preparation

`prepare_sucrose_spectra.py` is the canonical preparation step for the selected
600, 800, 900, and 1100 MHz sucrose spectra in `data/sucrose`.

It performs the same preparation for every field:

1. reads the processed Bruker `1r` and its `procs` metadata;
2. reconstructs the processed ppm axis and applies `NC_proc` scaling;
3. detects DSS within -0.20 to 0.20 ppm and shifts DSS to exactly 0 ppm;
4. defines the 3.00 to 5.80 ppm sucrose fitting region;
5. excludes water at 4.65 to 4.90 ppm;
6. excludes the unexplained feature at 5.15 to 5.30 ppm.

The 1100 MHz input uses `pdata/999`, the preferred 8x zero-filled
reprocessing. The other fields use `pdata/1`.

## Run

From the repository root:

```bash
python3 src/sucrose/prepare_sucrose_spectra.py
```

The script has no third-party Python dependencies. Outputs are written to
`outputs/sucrose/prepared`:

- `*_full.csv` contains the complete DSS-referenced spectrum, a baseline-
  corrected intensity, and the inclusion/exclusion label for every point.
- `*_fit.csv` contains only points retained for fitting.
- `preparation_summary.csv` records the Bruker metadata, detected DSS shift,
  anomeric peak position, and point counts for every field.

Run a subset with, for example:

```bash
python3 src/sucrose/prepare_sucrose_spectra.py --dataset 800
```

## Simulation and fitting source

- `simulate_sucrose_spinach_fft.m` is the shared Spinach forward model. It
  supports direct and noesypr1d acquisition modes and exact glucose/fructose
  seven-spin block splitting.
- `noesypr1d_acquire.m` is the repository-local noesypr1d pulse-sequence
  callback.
- `run_sucrose_field.m` configures each field using repository-relative paths.
- `capture_sucrose_fid.m` runs one configured field without plotting and saves
  the raw simulated complex FID before apodization and Fourier transformation.
- `bayes_astro/stock_transfer_test.py` fits the clean 600/900 MHz 5 mm-tube
  pair, then tests transfer to 800 MHz (same stock, different tube) and 1100
  MHz (different stock, same concentration).
- `bayes_astro/sucrose_sim.py` is an independent exact Python simulator with
  self-tests.

Run the Python physics self-tests before MATLAB:

```bash
python3 src/sucrose/bayes_astro/sucrose_sim.py
```

To capture and save a simulated sucrose FID from MATLAB, set
`SPINACH_ROOT`, change to the repository root, and run:

```matlab
addpath('/Users/cece/Desktop/edison_lab/final_repo/src/sucrose')
capture_sucrose_fid('600')
```

Use `'800'`, `'900'`, or `'1100'` for the other fields. The raw complex FID
is saved as both `.mat` and `.csv` under
`outputs/sucrose/<field>MHz/`. The CSV columns are time, real FID, imaginary
FID, and magnitude. A two-panel PNG of the full FID and its early-time
oscillations is saved alongside them. This is a simulated FID; the processed
Bruker `1r` file is the experimental frequency-domain spectrum.

Run the stock/tube transfer test after preparation:

```bash
python3 src/sucrose/bayes_astro/stock_transfer_test.py
```

It fits 600/900 MHz, then scores 800 and 1100 with nuisance parameters only.
Results are written to `outputs/sucrose/stock_transfer_600_900_to_800_1100.csv`.

After a candidate matrix passes the Python transfer gate, validate it with the
same Spinach/NOESY workflow without changing the configured GISSMO seed:

```bash
./run_spinach_field.sh sucrose --matrix \
  outputs/sucrose/sucrose_matrix_fit_600_900.txt
```

Candidate results are written under `outputs/sucrose/*MHz_spinach_candidate/`
and to `spinach_multifield_summary_candidate.csv`.

Compare the candidate with the original GISSMO Spinach runs:

```bash
conda run -n sucrose_project \
  python src/common/compare_spinach_matrices.py --molecule sucrose
```

The field-by-field changes are written to
`outputs/sucrose/spinach_candidate_vs_gissmo.csv`.

The automatic quality gate is implemented in
`src/common/multifield_quality_gate.py`. It reports:

- `PASS`: training and held-out fields meet the configured correlation, RMSE,
  transfer-improvement, and calibration checks;
- `REVIEW`: borderline data that should be inspected by a mentor;
- `REDO`: clearly failed transfer or training quality, suggesting reprocessing
  or reacquisition before changing the spin matrix.

Thresholds are in the molecule configuration under `quality_gate`, so the
same decision logic can be reused for another carbohydrate with a new matrix,
field manifest, and appropriate validation thresholds.

## Selecting a seed matrix

The project endpoint is a structure-specific GISSMO matrix that can be shared
with the carbohydrate-NMR community. GISSMO is preferred when an exact matrix
already exists. For a new carbohydrate, Bubb-style structural reporters and
the measured spectra are used to make a provisional seed:

```bash
python3 src/common/select_seed_matrix.py --molecule sucrose --seed-source auto --yes
```

Use `--interactive` when teaching and you want the student to choose among
available sources. The selector writes `outputs/<molecule>/seed_selection.json`.
If no spectra have been entered yet, the selector reports
`NEEDS_SPECTRAL_ASSIGNMENT`; that means “collect the observations needed to
build the seed,” not “abandon the molecule.” Build the seed with:

```bash
python3 src/common/build_provisional_seed.py --molecule new_sugar \
  --observations data/new_sugar/spectral_observations.json
```

That writes a `PROVISIONAL_SEED` matrix and provenance manifest. It may be used
for fitting and identifiability tests, but it cannot be labeled publishable
until the shared multi-field fit, held-out transfer, and uncertainty checks
pass. The final deposited artifact is the refined matrix plus its atom order,
source spectra, fit report, and provenance—not an untracked rule-based guess.

To import a BMRB entry first, use its metabolite entry ID (for example,
`bmse000119`):

```bash
python3 src/common/query_bmrb_entry.py --molecule sucrose --entry bmse000119
```

This uses BMRB's NMR-STAR REST endpoint, saves the raw entry, extracts the
assigned proton shifts, and writes `spectral_observations.json` plus a
`provenance.json`. It also checks the linked GISSMO entry. If a GISSMO project
exists, it downloads the project XML/archive, extracts the displayed matrix,
and writes a `READY` seed manifest pointing to that matrix. If no GISSMO
matrix exists, add the Bubb-guided connectivity and measured couplings to the
observation file before building the provisional matrix. BMRB is a shift
source; GISSMO is the preferred structure-specific matrix when available.

The current matrix-validation grouping is explicit: 600/900 MHz are the clean
same-stock, same-tube training pair; 800 MHz is the same stock in a different
tube; 1100 MHz is a different stock at the same concentration. The transfer
test keeps the latter two out of the physical-matrix training step while still
testing whether the shared matrix predicts them.

For a new carbohydrate with one to five fields, record `sample_id`, `tube_id`,
concentration, and (when known) a role for each dataset, then preview the
metadata-aware split:

```bash
python3 src/common/plan_multifield_split.py --molecule sucrose
```

The planner holds out entire sample/tube groups when possible, rather than
splitting by field number. One field gives `NO_MATRIX_VALIDATION`; two fields
give a provisional fit plus one holdout; three or more fields enable repeated
held-out validation. This prevents students from claiming transferability
when the training and validation spectra came from the same physical tube.

Students can import a new processed Bruker experiment without manually
copying files or editing JSON:

If a folder contains several consistently named experiments such as
`glucose_500MHz_..._exp3`, use the batch importer instead:

```bash
python3 src/common/import_bruker_bundle.py \
  --molecule glucose \
  --source-root /Users/cece/Desktop/glucose_BMRB_multifield \
  --auto-metadata
```

It detects all matching Bruker folders and imports them in one pass. The
auto mode reads concentration from the bundle README, reads pulse sequence
from `acqus`, assigns training/validation only among comparable conditions,
and leaves mismatched conditions `unassigned`. Students should review the
generated config before fitting.

For a completely new molecule, initialize its configuration first:

```bash
python3 src/common/query_bmrb_entry.py --molecule glucose --entry bmse000015
python3 src/common/init_carbohydrate.py --molecule glucose
```

The initializer uses the downloaded GISSMO seed when one is available.

```bash
python3 src/common/import_bruker_dataset.py \
  --molecule sucrose \
  --source /path/to/experiment \
  --field-mhz 1000 \
  --experiment 12 \
  --procno 1 \
  --sample-id stock_A \
  --tube-id 5mm \
  --concentration-mM 100 \
  --role validation
```

The importer requires `acqus`, `pdata/<procno>/1r`, and
`pdata/<procno>/procs`, copies them into `data/sucrose/1000_MHz/12`, and
updates `sucrose_config.json`. It refuses duplicate field keys unless
`--replace` is supplied.

## Workflow configuration

`data/sucrose/sucrose_config.json` is the single source of truth for settings
that describe the molecule and analysis rather than the spectrometer: the
matrix and atom order, independent spin blocks, fit/exclusion windows, pulse
sequence defaults, Spinach basis, and independent-fit bounds. MATLAB loads it
through `src/common/load_carbohydrate_config.m`; Python loads it through
`src/common/carbohydrate_config.py`. A new carbohydrate can copy this file,
replace the matrix/blocks/windows, and reuse the same preparation and fitting
code without editing field numbers into a script.

Acquisition values are intentionally not in this JSON. They are measured from
each Bruker `acqus`, `procs`, and `1r` file and regenerated into the preparation
summary. This keeps a sample-specific DSS correction (including the 800 MHz
offset) from becoming a magic number in the physics code.

## Metadata-driven acquisition setup

Do not hand-copy Bruker acquisition numbers into a driver. The source of each
number is:

- `acqus`: spectrometer frequency (`SFO1`), carrier (`O1`), sweep width
  (`SW_h`), acquired `TD`, pulse program, and scan count;
- `pdata/<procno>/procs`: processed size (`SI`), processed frequency (`SF`),
  processed sweep (`SW_p`), raw processed offset (`OFFSET`), byte order, data
  type, and `NC_proc` scaling;
- processed `1r`: the measured DSS position, used to calculate the
  DSS-referenced offset.

The reusable extractor is:

```bash
python3 src/common/bruker_metadata.py \
  --dataset data/sucrose/600_MHz/5 \
  --dataset data/sucrose/800_MHz/2 \
  --dataset data/sucrose/900_MHz/6 \
  --dataset data/sucrose/1100_MHz/7 \
  --procno 1 --procno 1 --procno 1 --procno 999 \
  --output outputs/sucrose/field_metadata.csv
```

`prepare_sucrose_spectra.py` writes the same metadata into
`outputs/sucrose/prepared/preparation_summary.csv`. The MATLAB field runner
reads that generated summary, including the DSS-referenced offset, so its
field-specific acquisition values are no longer embedded in the MATLAB code.
For another carbohydrate, point the extractor at that carbohydrate's Bruker
experiment directories and use the resulting metadata table in the analogous
runner.

## Optional nuisance diagnostic

After producing the fixed 600/900 MHz candidate matrix, test whether a small
additive baseline is limiting agreement without changing the matrix:

```bash
python3 src/sucrose/bayes_astro/nuisance_diagnostic.py \
  --candidate-matrix outputs/sucrose/sucrose_matrix_fit_600_900.txt
```

The script compares no-baseline, constant-baseline, and weak quadratic-
baseline fits for all four fields. It writes audit traces, metrics, and a plot
under `outputs/sucrose/nuisance_test/`. This is diagnostic only; accept a
correction only if the blind 800/1100 MHz results improve.
