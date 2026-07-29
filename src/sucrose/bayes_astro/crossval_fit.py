"""
crossval_fit.py -- Stage 2: leave-one-field-out cross-validation of the joint
matrix refinement.

WHY: joint_fit.py improved all three fields' agreement, but a joint fit forces
ONE shared matrix and therefore CANNOT distinguish refinement from overfitting --
degenerate crowded-region params simply absorb model mismatch. (Proof: it
inflated the isolated-CH2 splitting to 0.012 ppm, a value we ALREADY showed is a
field artifact via the per-field three-field test.)

THE HONEST TEST: does a refinement learned on some fields PREDICT a field it
never saw? For each held-out field f:
    - refine the full matrix on the OTHER two fields (blind to f)
    - then, for f, fit ONLY its instrument nuisances (lb, offset) -- those are
      legit per-field acquisition params, not the physics under test -- and score.
Compare three matrices, each scored on f with f's own best lb/offset:
    r_GISSMO(f)  -- deposited matrix, no shift/coupling fitting   (baseline)
    r_CV(f)      -- matrix refined on the OTHER two, PREDICTING f  (the real test)
    r_train(f)   -- matrix refined on ALL three, incl f (in-sample upper bound)

READING:
    r_CV > r_GISSMO on ALL held-out fields  => the multi-field refinement
        genuinely transfers to unseen fields -> real, keep it.
    r_CV ~ or < r_GISSMO                     => the joint gains are overfitting;
        the GISSMO deposit is the better community matrix.
The CH2-split column shows each BLIND fit's isolated-CH2 splitting: a REAL value
would be ~consistent row-to-row; scatter (or all >> GISSMO's 0.003) => artifact.

Usage: python3 crossval_fit.py     (~4 short least-squares fits, a few minutes)
Requires: numpy, scipy, and joint_fit.py in the same dir.
"""
import numpy as np
from scipy.optimize import least_squares
import joint_fit as jf


LB0 = dict(zip(jf.NAMES, jf.MODEL_CONFIG["initial_linewidth_hz"]))
NCOUP = len(jf.CIDX)
NSPIN = jf.N_SPINS
SHIFT_BOUND = float(jf.MODEL_CONFIG["shift_bound_ppm"])
COUPLING_BOUND = float(jf.MODEL_CONFIG["coupling_bound_hz"])
LB_BOUNDS = tuple(jf.MODEL_CONFIG["linewidth_bounds_hz"])
OFFSET_BOUND = float(jf.MODEL_CONFIG["offset_bound_ppm"])


def _build_M(theta):
    M = np.diag(theta[:NSPIN]).astype(float)
    for (i, j), c in zip(jf.CIDX, theta[NSPIN:NSPIN + NCOUP]):
        M[i, j] = c
    return M


def fit_on(train_names):
    """Refine the full matrix (+ ONLY the training fields' nuisances).

    Crucial: we do NOT carry the held-out field's lb/offset in the parameter
    vector. If we did, those params would have zero gradient (they touch no
    residual) -> a zero-norm Jacobian column -> x_scale='jac' divides by zero
    and trf stalls. Parametrizing exactly the free params keeps it well-scaled
    and fast (same as the all-three joint fit)."""
    coup0 = np.array([jf.J0[i, j] for (i, j) in jf.CIDX])
    nuis0 = np.ravel([[LB0[n], 0.0] for n in train_names])
    theta0 = np.concatenate([jf.SHIFTS0, coup0, nuis0])
    lo = np.concatenate([jf.SHIFTS0 - SHIFT_BOUND, coup0 - COUPLING_BOUND,
                         np.ravel([[LB_BOUNDS[0], -OFFSET_BOUND] for _ in train_names])])
    hi = np.concatenate([jf.SHIFTS0 + SHIFT_BOUND, coup0 + COUPLING_BOUND,
                         np.ravel([[LB_BOUNDS[1], OFFSET_BOUND] for _ in train_names])])

    def resid(theta):
        M = _build_M(theta)
        rest = theta[NSPIN + NCOUP:]
        out = []
        for k, name in enumerate(train_names):
            d = jf.DATA[name]
            lb, off = rest[2 * k], rest[2 * k + 1]
            out.append((jf.field_model(M, d, lb, off) - d["y"]) / d["sig"])
        return np.concatenate(out)

    # Optimizer notes (a config sweep settled this): x_scale='jac' with the exact
    # trust-region solver made trf CRAWL/stall at cost ~2.1e4 (ill-conditioned,
    # degenerate 2-field problem); loosening tol instead caused erratic premature
    # xtol stops (some folds quit in ~15 evals at very different quality).
    # tr_solver='lsmr' (iterative, tolerant of the rank-deficient Jacobian)
    # converges cleanly via ftol at consistent quality -- keep tight default tols.
    sol = least_squares(resid, theta0, bounds=(lo, hi), method="trf",
                        tr_solver="lsmr", max_nfev=int(jf.MODEL_CONFIG.get("crossval_max_nfev", 1000)), verbose=1)
    return _build_M(sol.x)             # return the refined matrix


def best_nuis_r(M, name):
    """Score matrix M on field `name`, giving that field its OWN best lb/offset."""
    d = jf.DATA[name]

    def resid(p):
        return (jf.field_model(M, d, p[0], p[1]) - d["y"]) / d["sig"]

    sol = least_squares(resid, [LB0[name], 0.0],
                        bounds=([LB_BOUNDS[0], -OFFSET_BOUND], [LB_BOUNDS[1], OFFSET_BOUND]),
                        max_nfev=int(jf.MODEL_CONFIG.get("nuisance_max_nfev", 800)))
    sim = jf.field_model(M, d, sol.x[0], sol.x[1])
    g = np.isfinite(sim) & np.isfinite(d["y"])
    return np.corrcoef(sim[g], d["y"][g])[0, 1], sol.x[0], sol.x[1]


def main():
    theta0, _ = jf.theta0_and_bounds()
    M_gissmo = jf.unpack(theta0)[0]

    print("\nfitting all-three (in-sample upper bound) ...")
    M_all = fit_on(jf.NAMES)

    print("\n===== leave-one-field-out cross-validation =====")
    print(" held-out |  r_GISSMO   r_CV(PREDICT)   r_train(in-sample) | CH2 split | verdict")
    print("          |             <-- real test -->                 | (GISSMO .0030)")
    n_gen = 0
    for held in jf.NAMES:
        train = [n for n in jf.NAMES if n != held]
        print(f"\n  ... refining on {train} (blind to {held}) ...")
        M_cv = fit_on(train)
        r_g, _, _ = best_nuis_r(M_gissmo, held)
        r_cv, _, _ = best_nuis_r(M_cv, held)
        r_tr, _, _ = best_nuis_r(M_all, held)
        split_cv = M_cv[13, 13] - M_cv[12, 12]
        if r_cv > r_g + 0.01:
            verdict, ok = "generalizes", True
        elif r_cv > r_g - 0.01:
            verdict, ok = "neutral", False
        else:
            verdict, ok = "OVERFIT", False
        n_gen += ok
        print(f"  {held:7s}  |  {r_g:.4f}     {r_cv:.4f}          {r_tr:.4f}       "
              f"| {split_cv:+.4f}  | {verdict}")

    print("\n" + "=" * 66)
    if n_gen == len(jf.NAMES):
        print("VERDICT: the joint refinement improves EVERY held-out field it never")
        print("saw -> it generalizes. The refined matrix is a real improvement;")
        print("promote survivors to Stage-3 MCMC for uncertainties.")
    else:
        print("VERDICT: the joint refinement does NOT reliably beat GISSMO on")
        print("held-out fields -> the improved in-sample r was largely overfitting")
        print("the crowded region. GISSMO's deposit stands as the community matrix;")
        print("1D data alone can't refine the aliphatic pileup. (This is itself the")
        print("publishable result: field-transferability, not fit quality, is the")
        print("criterion -- and it says 'don't touch the deposit'.)")
    print("Watch the CH2-split column: real => steady row-to-row; scattered or all")
    print(">> 0.003 => the isolated CH2 is absorbing mismatch (known artifact).")


if __name__ == "__main__":
    main()
