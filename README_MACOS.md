# macOS setup

This guide is for students using macOS. Run the commands from Terminal.

## 1. Download the repository

Open the repository in a browser:
<https://github.com/edisonomics/carbohydrate-nmr-spin-matrices>. Select the
green **Code** button, choose **Download ZIP**, and double-click the ZIP in
your Downloads folder. Move the extracted folder to your Desktop if desired,
then open Terminal and run (adjust the path if you chose another location):

```bash
cd ~/Desktop/carbohydrate-nmr-spin-matrices-main
```

macOS Terminal normally uses **zsh**. That is fine for this project: run the
commands in zsh, and execute the included Bash scripts directly. You do not
need to switch shells or type `bash` before each command.

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

MATLAB and the Spinach-required toolboxes must be installed separately. For
Spinach 2.10.1, the tested setup requires **MATLAB R2024a or later** (the full
desktop product; MATLAB Runtime is not sufficient) and:

- **Parallel Computing Toolbox**
- **Deep Learning Toolbox**
- **Reinforcement Learning Toolbox**
- **Optimization Toolbox**
- **Statistics and Machine Learning Toolbox**
- **Mapping Toolbox**

The repository installer downloads the tested Spinach 2.10.1 release into the
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
on `PATH`, set `MATLAB_BIN` to its executable before running the launcher. For
example, for a standard MATLAB R2024a application installation:

```bash
export MATLAB_BIN="/Applications/MATLAB_R2024a.app/bin/matlab"
./run_spinach_field.sh sucrose
```

Change `R2024a` to the release installed on your Mac. You can check the path
with `ls /Applications/MATLAB_*.app/bin/matlab`.
