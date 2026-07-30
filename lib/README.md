# Spinach dependency

Spinach is a third-party MATLAB library and is intentionally not distributed
with this repository. Install Spinach separately. When it is installed in
the repository's standard location, the project scripts set `SPINACH_ROOT`
automatically.

Recommended repository installation:

```bash
./scripts/install_spinach.sh
source ./setup_paths.sh
```

`setup_paths.sh` sets:

```text
SPINACH_ROOT=<repository>/lib/Spinach-2.10.1
```

Only use a manual override when Spinach is installed somewhere else:

```bash
export SPINACH_ROOT="$HOME/Spinach-2.10.1"
```
