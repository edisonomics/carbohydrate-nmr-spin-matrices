# Edison Lab spin-matrix repository

Load the repository paths before running MATLAB workflows:

```bash
source ./setup_paths.sh
```

The script exports `SPINACH_ROOT`, `EDISON_DATA_ROOT`,
`EDISON_OUTPUT_ROOT`, `ALANINE_SRC`, and `SUCROSE_SRC` relative to this
repository, so the setup works regardless of the directory from which it is
sourced.

## Dependencies

Install these before running the workflow:

- **MATLAB** with the toolboxes required by Spinach (including Parallel
  Computing and Optimization Toolbox).
- **Spinach 2.10.1**, downloaded separately from the official
  [Spinach releases](https://github.com/IlyaKuprov/Spinach/releases).
- **Miniconda**, which installs Python and the pinned packages in
  [`environment.sucrose.yml`](environment.sucrose.yml).

### Install Miniconda (one time)

If Terminal says `conda: command not found`, install Miniconda from the
official [Miniconda download page](https://docs.conda.io/miniconda.html).
On macOS, open Terminal and run `uname -m`: choose **Apple Silicon/arm64** if
it prints `arm64`, or **Intel/x86_64** if it prints `x86_64`. The graphical
installer is the simplest choice. Allow the installer to initialize your
shell, then close and reopen Terminal.

Confirm that it worked:

```bash
conda --version
```

From the repository root, create this project's environment:

```bash
conda env create -f environment.sucrose.yml
conda activate sucrose
```

You only need to create the environment once. In later sessions, use
`conda activate sucrose` before running the Python commands.

### Windows note

Windows students should use the **Anaconda Prompt** installed with Miniconda.
The Python workflow is the same; use Windows paths and replace `/` with `\\`
if needed:

```text
cd C:\path\to\carbohydrate-nmr-spin-matrices
conda activate sucrose
python src\sucrose\prepare_sucrose_spectra.py
```

The `.sh` launcher and `source ./setup_paths.sh` require a Unix-like shell.
For Spinach, either run them through **WSL2** or **Git Bash**, or set the
Spinach path in the Anaconda Prompt before launching MATLAB:

```bat
set SPINACH_ROOT=C:\path\to\carbohydrate-nmr-spin-matrices\lib\Spinach-2.10.1
```

The path must contain `kernel`, `etc`, `experiments`, and `interfaces`.

Spinach is a third-party MATLAB library and is not included in this Git
repository. Download and extract the `Spinach-2.10.1` release **inside this
repository's `lib/` folder**. The expected layout is:

```text
carbohydrate-nmr-spin-matrices/
└── lib/
    └── Spinach-2.10.1/
        ├── kernel/
        ├── etc/
        ├── experiments/
        └── interfaces/
```

If `lib/Spinach-2.10.1/kernel` exists, the installation is in the right
place. If the archive creates an extra nested folder, move the inner
`Spinach-2.10.1` directory up one level. Alternatively, override
`SPINACH_ROOT` after loading the repository paths for an installation
elsewhere:

```bash
source ./setup_paths.sh
export SPINACH_ROOT="$HOME/Spinach-2.10.1"
```

The `lib/` directory is intentionally ignored by Git, so the Spinach source
and local experimental data are never committed. See the
[Spinach installation guide](https://spindynamics.org/wiki/index.php?title=Installation)
for MATLAB path and version requirements.

## Prepare the official sucrose spectra

The reproducible Bruker `1r` loading, DSS referencing, and fit-window masking
workflow for the 600, 800, 900, and 1100 MHz sucrose data is:

```bash
python3 src/sucrose/prepare_sucrose_spectra.py
```

Prepared full and fit-only spectra are written to
`outputs/sucrose/prepared`. See `src/sucrose/README.md` for details.

## Run one Spinach field

Students do not need to set MATLAB paths or choose a field manually. From the
repository root, run all prepared fields with:

```bash
./run_spinach_field.sh glucose
```

The launcher finds the bundled Spinach library and MATLAB, then runs the
generic field wrapper once per prepared dataset. Results are written under
`outputs/glucose/*MHz_spinach`.

## Record carbohydrate chemistry evidence

The seed workflow combines molecule-specific BMRB/GISSMO artifacts with
machine-readable guidance from Bubb (2003) and companion assignment/coupling
references. It flags a flattened single-spin model when a reducing sugar
should be represented as separate anomer components:

```bash
python3 src/common/assess_carbohydrate_evidence.py --molecule glucose
```

The report is written to
`outputs/<molecule>/chemistry_evidence.json`. `PASS` means the configured
model is chemically consistent with the selected profile; `REVIEW` means a
student should resolve the warning before treating the fit as publishable.

The BMRB importer also preserves deposited numeric peak coordinates and atom
labels in `data/<molecule>/bmrb/<entry>/spectral_observations.json` under
`peak_data`. It extracts any available 1-D or 2-D `_Peak_char` tables (for
example, HSQC cross-peaks) but never turns peak positions into scalar
couplings automatically. An experiment can therefore be documented even
when BMRB does not provide its numeric peak table.

For an explicit, non-destructive HSQC diagnostic on a reducing-sugar mixture:

```bash
python3 src/common/refine_hsqc_constraints.py --molecule mannose
```

This writes a candidate and a before/after multifield comparison under
`outputs/<molecule>/hsqc_refinement`; it never replaces the configured matrix.

### Independent J-measurement stage

Chemical shifts and 2-D cross-peaks establish where the spins are and how
they connect; they do not, by themselves, provide signed scalar couplings.
Before refining a new carbohydrate matrix, use resolved 1-D multiplet
splittings (or a J-resolved 2-D experiment) as an independent J check.  The
screening tool is intentionally conservative: it reports candidate spacings,
records the raw-FID provenance, and leaves the matrix unchanged.

For example, after importing a provisional xylose seed and its 1-D spectra:

```bash
conda run -n sucrose_project \
  python src/common/measure_1d_j_candidates.py \
  --molecule xylose --fields 600 900 1100
```

The report is written to
`outputs/<molecule>/j_measurements/j_candidate_measurements.json` and is
always marked `REVIEW` until a person assigns each multiplet and confirms the
same coupling across fields.  A spacing can be used to update a provisional
matrix only after that review; it must never overwrite a verified GISSMO
matrix automatically.  If no line is resolved, retain the Bubb/BMRB value as
a documented prior and record that an independent J measurement was
unavailable.  `assess_carbohydrate_evidence.py` now includes this stage in
`chemistry_evidence.json` as `j_measurement_evidence` (`NOT_RUN` or `REVIEW`).

For troubleshooting one field only, an advanced invocation is available:

```bash
./run_spinach_field.sh glucose 700
```
