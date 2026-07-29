# Alanine spin-matrix code

This directory contains the reusable alanine AX₃ Spinach workflow. It uses the
lab data in `final_repo/data/alanine` and writes generated outputs to
`final_repo/outputs/alanine`.

## Run a field-transfer fit

Set the Spinach installation once:

```bash
export SPINACH_ROOT=/path/to/Spinach-2.10.1
```

Then run from MATLAB:

```matlab
run('src/alanine/run_alanine_600MHz.m')
run('src/alanine/run_alanine_800MHz.m')
run('src/alanine/run_alanine_900MHz.m')
```

The primary lab-sample comparison is 600/800/900 MHz. The downloaded 700 MHz
BMRB spectrum is intentionally excluded.

The drivers call `simulate_alanine_spinach_fft.m`, which builds the four-spin AX₃
model and fits linewidth plus a constant ppm calibration offset against the
processed Bruker `1r` spectrum.

For the full Bruker-style noesypr1d validation fits, run:

```matlab
run('src/alanine/fit_alanine_noesypr1d_experiment.m')
run('src/alanine/fit_alanine_noesypr1d_800MHz.m')
run('src/alanine/fit_alanine_noesypr1d_900MHz.m')
```

These keep the AX₃ matrix fixed and write field-specific results under
`final_repo/outputs/alanine/{600MHz,800MHz,900MHz}`.

To regenerate the Georgia-style 600/900 MHz comparison figures after running
the field drivers:

```bash
conda run -n sucrose_project python src/alanine/make_alanine_georgia_plots.py
```

The figures are written to `final_repo/outputs/alanine/georgia_plots/`. Each
field has a full alanine-window plot and a two-panel CH₃/Hα zoom. The same
command also writes combined 600/900 MHz poster figures named
`alanine_600_900MHz_georgia_full_2panel.*` and
`alanine_600_900MHz_georgia_zooms_2panel.*`. The traces
are the experimental spectrum, the direct Spinach FFT, and the analytical
first-order AX₃ theory using the shared 3.7680/1.4655 ppm shifts and 7.234 Hz
coupling.

For the actual Bruker-style noesypr1d figures, rerun the 600/900 MHz scripts
and then format the three-way overlays:

```bash
conda run -n sucrose_project python src/alanine/make_alanine_noesypr1d_georgia_plots.py
```

Those outputs are written to `final_repo/outputs/alanine/georgia_noesypr1d_plots/`
and include experiment, noesypr1d Spinach fit, and AX₃ theory.

To diagnose the 800 MHz acquisition mismatch, compare all three imported
800 MHz acquisitions with one shared Spinach simulation:

```matlab
run('src/alanine/fit_alanine_noesypr1d_800MHz_all.m')
```

This writes `outputs/alanine/800MHz_diagnostic/alanine_800MHz_acquisition_comparison.csv`
and a comparison figure. It reads each acquisition's Bruker `procs` metadata,
including its own `NC_proc` scaling, so acquisition 1/8 are not silently
processed as acquisition 11. The diagnostic estimates a starting calibration
offset from the observed methyl and Halpha anchors, rather than assuming that
every field is already close to zero ppm offset. Choose an acquisition with a
physical linewidth (roughly 1–3 Hz) and low RMSE before using it as the 800 MHz
validation trace.
