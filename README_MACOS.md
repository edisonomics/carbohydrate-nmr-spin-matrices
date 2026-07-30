# macOS setup

This guide is for students using macOS. Run the commands from Terminal.

## 1. Clone the repository

```bash
git clone https://github.com/edisonomics/carbohydrate-nmr-spin-matrices.git
cd carbohydrate-nmr-spin-matrices
```

## 2. Install Miniconda

Download the graphical installer from the official
[Miniconda page](https://docs.conda.io/miniconda.html). In Terminal, run
`uname -m` to identify the installer:

- `arm64` means Apple Silicon (M1/M2/M3/M4).
- `x86_64` means Intel.

Allow the installer to initialize your shell, then close and reopen Terminal.
Verify the installation:

```bash
conda --version
```

Create and activate the project environment:

```bash
conda env create -f environment.sucrose.yml
conda activate sucrose
```

## 3. Install Spinach

MATLAB and the Spinach-required toolboxes must be installed separately. The
repository installer downloads the tested Spinach 2.10.1 release into the
ignored `lib/` directory:

```bash
./scripts/install_spinach.sh
source ./setup_paths.sh
```

The expected check is:

```bash
test -d lib/Spinach-2.10.1/kernel && echo "Spinach is ready"
```

## 4. Run the Python workflow

```bash
conda activate sucrose
python3 src/sucrose/prepare_sucrose_spectra.py
python3 -m unittest discover -s tests -v
```

Place local Bruker data under `data/<molecule>/`; experimental data and
generated `outputs/` are intentionally ignored by Git.

## 5. Run Spinach

After MATLAB is installed and `source ./setup_paths.sh` has been run:

```bash
./run_spinach_field.sh sucrose
```

The launcher uses the MATLAB executable found on the system. If MATLAB is not
on `PATH`, set `MATLAB_BIN` to its executable before running the launcher.
