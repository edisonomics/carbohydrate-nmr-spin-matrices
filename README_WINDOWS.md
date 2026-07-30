# Windows setup through WSL2

Use WSL2 (Windows Subsystem for Linux 2) for this project. It provides the
Linux shell expected by the repository's Bash scripts and avoids mixing
Windows and Linux Python environments.

## 1. Install WSL2

Open **PowerShell as Administrator** and run:

```powershell
wsl --install -d Ubuntu
```

Restart Windows when prompted. Launch **Ubuntu** from the Start menu and
create the Linux username and password it requests. Microsoft maintains the
[WSL installation guide](https://learn.microsoft.com/windows/wsl/install).

Install the small command-line utilities used by the setup scripts:

```bash
sudo apt update
sudo apt install -y git curl unzip
```

## 2. Install Miniconda inside Ubuntu

Do not use the Windows Miniconda installation for the WSL environment. In the
Ubuntu window, follow the Linux instructions on the official
[Miniconda page](https://docs.conda.io/miniconda.html), then close and reopen
Ubuntu. Confirm that Conda is available:

```bash
conda --version
```

## 3. Clone the repository inside WSL2

Keep the working copy in the Linux home directory for better performance:

```bash
cd ~
git clone https://github.com/edisonomics/carbohydrate-nmr-spin-matrices.git
cd carbohydrate-nmr-spin-matrices
```

Create the project environment:

```bash
conda env create -f environment.sucrose.yml
conda activate sucrose
```

## 4. Install Spinach inside WSL2

```bash
./scripts/install_spinach.sh
source ./setup_paths.sh
test -d lib/Spinach-2.10.1/kernel && echo "Spinach is ready"
```

## 5. Run the Python workflow

```bash
conda activate sucrose
python3 src/sucrose/prepare_sucrose_spectra.py
python3 -m unittest discover -s tests -v
```

If your data is on the Windows desktop, it is visible under `/mnt/c`, for
example `/mnt/c/Users/YourName/Desktop/my_nmr_data`. Copying it into the WSL
home directory is usually faster:

```bash
cp -R /mnt/c/Users/YourName/Desktop/my_nmr_data data/my_molecule
```

## 6. MATLAB and Spinach on Windows

WSL2 does not install MATLAB. The Spinach stage requires a licensed MATLAB
installation. You may either install MATLAB for Linux inside WSL2, or use the
Windows MATLAB application separately for the `.m` workflows. If MATLAB is
installed only on Windows, complete the Python workflow in WSL2 and run the
MATLAB scripts from the Windows MATLAB interface using the repository path.

The Bash launcher is intended for a MATLAB executable available inside WSL2:

```bash
./run_spinach_field.sh sucrose
```
