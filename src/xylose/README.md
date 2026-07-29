# Xylose multifield fitting

The xylose data are a same-tube 5 mm multifield series. The configured split
uses 600 and 900 MHz as training fields and holds 1100 MHz out for validation.

Run the shared alpha/beta physical fit from the repository root:

```bash
conda run -n sucrose_project python src/xylose/joint_mixture_fit.py
```

The fit shares 12 chemical shifts, 12 nonzero J couplings, and the alpha
fraction across fields. Alpha/beta linewidths, calibration offset, scale, and
baseline are field-specific nuisance parameters. The 1100 MHz spectrum never
changes the physical matrices; it is used only for held-out transfer scoring.

Before calling the provisional matrix a publishable new matrix, screen the
resolved 1-D multiplets for independent J evidence:

```bash
conda run -n sucrose_project \
  python src/common/measure_1d_j_candidates.py \
  --molecule xylose --fields 600 900 1100
```

This is a review aid, not an automatic refiner.  The candidate spacings are
stored under `outputs/xylose/j_measurements`; they must be assigned manually
and checked for cross-field agreement before entering the matrix.  COSY and
TOCSY support connectivity/assignment, while HSQC supports one-bond H-C
assignment; neither should be reported as a numeric J measurement unless the
experiment is explicitly J-resolved.

Run all three two-field training folds with one blind field each:

```bash
conda run -n sucrose_project python src/xylose/crossval_mixture_fit.py
```

Fold-specific matrices and blind-prediction metrics are written under
`outputs/xylose/cross_validation`.
