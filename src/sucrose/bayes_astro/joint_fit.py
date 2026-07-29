"""
joint_fit.py -- Stage 1 of the full multi-field matrix analysis: a JOINT
maximum-a-posteriori (least-squares) fit of the ENTIRE sucrose matrix (all 14
shifts + all nonzero J-couplings) to the configured training fields (currently
600 and 900 MHz), each field with its own linewidth + calibration offset.

Why MAP (optimizer) not MCMC here: ~27 physical + 6 nuisance params is far too
high-dim for emcee to converge in reasonable time, and it would hit label
degeneracies for many near-equivalent protons. A local optimizer started at
GISSMO refines in place (no mode-swapping), runs in minutes, and directly shows
WHICH parameters want to move and by how much. Reproducibility across fields
(not any single fit) is then the arbiter of what's real; MCMC on the survivors
comes later for uncertainty.

Usage:
  python3 joint_fit.py --check     # evaluate at GISSMO, report per-field r/RMSE
                                   # + save overlays. VERIFY THIS FIRST.
  python3 joint_fit.py             # run the joint MAP fit
Requires: numpy, scipy, matplotlib.
"""
import os
import argparse
import csv
import numpy as np
from scipy.optimize import least_squares

import sucrose_sim as ss

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "common"))
from carbohydrate_config import load_config
from seed_manifest import resolve_seed_matrix

HERE = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATA_ROOT = os.path.join(REPO_ROOT, "data", "sucrose")
CONFIG = load_config(__import__("pathlib").Path(REPO_ROOT), "sucrose")
PROCESSING = CONFIG["processing"]
MODEL_CONFIG = CONFIG["independent_model"]
SUMMARY_FILE = os.path.join(REPO_ROOT, "outputs", "sucrose", "prepared",
                            "preparation_summary.csv")
PLAN_FILE = os.path.join(REPO_ROOT, "outputs", "sucrose", "multifield_plan.json")


def load_preparation_summary(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Missing {path}; run src/sucrose/prepare_sucrose_spectra.py first."
        )
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {str(row["field_mhz"]): row for row in rows}


SUMMARY = load_preparation_summary(SUMMARY_FILE)

def field_from_metadata(field_key):
    row = SUMMARY[str(field_key)]
    return dict(
        SFO1=float(row["sfo1_mhz"]), O1=float(row["o1_hz"]),
        exp=os.path.join(DATA_ROOT, row["relative_dir"], "pdata", row["procno"], "1r"),
        SI=int(row["points"]), SF=float(row["sf_mhz"]), SW_p=float(row["sw_hz"]),
        OFFSET=float(row["offset_dss_ppm"]), NC=int(row["nc_proc"]),
    )


# Fields are selected by the metadata-aware plan when available. Their
# acquisition numbers come from the generated Bruker metadata summary, never
# from this fit script. The `n` suffix preserves historical report labels.
if os.path.isfile(PLAN_FILE):
    with open(PLAN_FILE, encoding="utf-8") as handle:
        _PLAN = __import__("json").load(handle)
    _TRAINING_KEYS = [str(key) for key in _PLAN.get("training_fields", [])]
else:
    _TRAINING_KEYS = [str(key) for key in MODEL_CONFIG["fields_for_joint_fit"]]
FIELD_KEYS = {
    (f"{key}n" if str(key) in {"600", "1100"} else str(key)): str(key)
    for key in _TRAINING_KEYS
}
FIELDS = {label: field_from_metadata(key) for label, key in FIELD_KEYS.items()}
FIT_LO, FIT_HI = PROCESSING["fit_region_ppm"]
WATER, ARTIFACT, ANOM = (tuple(PROCESSING["water_region_ppm"]),
                         tuple(PROCESSING["artifact_region_ppm"]),
                         tuple(PROCESSING["anomeric_region_ppm"]))
STRIDE = int(MODEL_CONFIG["fit_stride"])

MATRIX0 = ss.load_gissmo_matrix(str(resolve_seed_matrix(__import__("pathlib").Path(REPO_ROOT), "sucrose", CONFIG)))
N_SPINS = MATRIX0.shape[0]
BLOCKS = [[int(index) - 1 for index in block] for block in CONFIG.get("blocks", [])]
SHIFTS0 = np.diag(MATRIX0).copy()
J0 = np.triu(MATRIX0, 1); J0 = J0 + J0.T
CIDX = [(i, j) for i in range(N_SPINS) for j in range(i + 1, N_SPINS)
        if J0[i, j] != 0.0]  # coupling positions
print(f"{len(SHIFTS0)} shifts, {len(CIDX)} nonzero couplings to fit.")


def load_field(fc):
    raw = np.fromfile(fc["exp"], dtype="<i4").astype(float) * (2.0 ** fc["NC"])
    n = raw.size
    ppm = fc["OFFSET"] - np.arange(n) * (fc["SW_p"] / fc["SF"]) / (n - 1)
    o = np.argsort(ppm); ppm, y = ppm[o], raw[o]
    y = y - np.median(y)
    fit = ((ppm >= FIT_LO) & (ppm <= FIT_HI)
           & ~((ppm >= WATER[0]) & (ppm <= WATER[1]))
           & ~((ppm >= ARTIFACT[0]) & (ppm <= ARTIFACT[1])))
    pf, yf = ppm[fit][::STRIDE], y[fit][::STRIDE]
    anom = (ppm >= ANOM[0]) & (ppm <= ANOM[1])
    yfn = yf / np.max(y[anom])
    amask = (pf >= ANOM[0]) & (pf <= ANOM[1])
    # Equal per-field weight: the dominant error is MODEL mismatch (~few % of
    # anomeric, similar across fields), not baseline noise (which differs ~7x
    # between fields and would let the lowest-noise field dominate). A common
    # sig weights the 3 fields equally; its exact value doesn't move the
    # least-squares minimum, only the reported chi-square scale.
    return dict(ppm=pf, y=yfn, anom=amask, sig=float(MODEL_CONFIG["noise_sigma"]),
                SFO1=fc["SFO1"], carrier=fc["O1"] / fc["SFO1"])


DATA = {name: load_field(fc) for name, fc in FIELDS.items()}
NAMES = list(FIELDS)


def unpack(theta):
    shifts = theta[:N_SPINS]
    coups = theta[N_SPINS:N_SPINS + len(CIDX)]
    rest = theta[N_SPINS + len(CIDX):]             # per-field [lb, offset] x nfields
    field_pars = {NAMES[k]: (rest[2 * k], rest[2 * k + 1]) for k in range(len(NAMES))}
    M = np.diag(shifts).astype(float)
    for (i, j), c in zip(CIDX, coups):
        M[i, j] = c
    return M, field_pars


def field_model(M, d, lb, offset):
    sppm, sint = ss.sucrose_sticks(M, d["SFO1"], d["carrier"], blocks=BLOCKS)
    sim = ss.lorentzian_spectrum(d["ppm"], sppm, sint, lb, d["SFO1"], offset)
    smax = np.max(sim[d["anom"]])
    if not np.isfinite(smax) or smax <= 0:
        return np.full_like(d["ppm"], 1e3)
    return sim / smax


def residuals(theta):
    M, fp = unpack(theta)
    res = []
    for name in NAMES:
        d = DATA[name]
        lb, off = fp[name]
        sim_n = field_model(M, d, lb, off)
        res.append((sim_n - d["y"]) / d["sig"])
    return np.concatenate(res)


def theta0_and_bounds():
    lb0 = list(MODEL_CONFIG["initial_linewidth_hz"])
    off0 = [0.0] * len(NAMES)
    theta0 = np.concatenate([SHIFTS0, [J0[i, j] for (i, j) in CIDX],
                             np.ravel(list(zip(lb0, off0)))])
    shift_bound = float(MODEL_CONFIG["shift_bound_ppm"])
    coupling_bound = float(MODEL_CONFIG["coupling_bound_hz"])
    lb_lo, lb_hi = MODEL_CONFIG["linewidth_bounds_hz"]
    off_bound = float(MODEL_CONFIG["offset_bound_ppm"])
    lo = np.concatenate([SHIFTS0 - shift_bound,
                         [J0[i, j] - coupling_bound for (i, j) in CIDX],
                         np.ravel(list(zip([lb_lo] * len(NAMES), [-off_bound] * len(NAMES))))])
    hi = np.concatenate([SHIFTS0 + shift_bound,
                         [J0[i, j] + coupling_bound for (i, j) in CIDX],
                         np.ravel(list(zip([lb_hi] * len(NAMES), [off_bound] * len(NAMES))))])
    return theta0, (lo, hi)


def report(theta, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    M, fp = unpack(theta)
    print(f"\n===== per-field agreement ({tag}) =====")
    fig, axes = plt.subplots(len(NAMES), 1, figsize=(14, 9), sharex=True)
    for ax, name in zip(axes, NAMES):
        d = DATA[name]; lb, off = fp[name]
        sim_n = field_model(M, d, lb, off)
        good = np.isfinite(sim_n) & np.isfinite(d["y"])
        r = np.corrcoef(sim_n[good], d["y"][good])[0, 1]
        rmse = np.sqrt(np.mean((sim_n[good] - d["y"][good]) ** 2))
        print(f"  {name:6s} r={r:.4f}  RMSE={rmse:.4f}  lb={lb:.2f} Hz  off={off:+.4f}")
        ax.plot(d["ppm"], d["y"], "k", lw=0.7, label=f"{name} exp")
        ax.plot(d["ppm"], sim_n, "r", lw=0.7, alpha=0.8, label="model")
        ax.legend(fontsize=8, loc="upper left"); ax.set_ylabel("norm")
    axes[-1].set_xlabel("1H shift (ppm)"); axes[0].invert_xaxis()
    fig.suptitle(f"Joint multi-field fit ({tag})")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, f"joint_fit_{tag}.png"), dpi=150)
    print(f"  wrote joint_fit_{tag}.png")
    return M


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="evaluate at GISSMO only")
    args = ap.parse_args()

    theta0, bounds = theta0_and_bounds()

    if args.check:
        report(theta0, "GISSMO")
        print("\nThis is the forward model at the deposited GISSMO matrix. Confirm the")
        print("per-field overlays look sane before running the fit (drop --check).")
        return

    print("\nrunning joint MAP fit (least_squares) over "
          f"{len(theta0)} params across {len(NAMES)} fields ...")
    # ~33 params -> each numerical-Jacobian step is ~34 residual evals, so allow
    # plenty of budget (each eval ~30-50 ms; a few thousand -> a few minutes).
    sol = least_squares(residuals, theta0, bounds=bounds, method="trf",
                        x_scale="jac", verbose=2, max_nfev=int(MODEL_CONFIG["max_nfev"]))
    M_fit = report(sol.x, "MAP")

    # ---- which parameters moved, and by how much ----
    print("\n===== shift deviations from GISSMO (Hz at 600 MHz) =====")
    sh = sol.x[:N_SPINS]
    for k in range(N_SPINS):
        d_ppm = sh[k] - SHIFTS0[k]
        reference_mhz = float(CONFIG["identifiability"]["reference_field_mhz"])
        hz = abs(d_ppm) * reference_mhz
        flag = "  <-- moved" if hz > float(CONFIG["reporting"]["shift_flag_hz"]) else ""
        print(f"  shift[{k:2d}] GISSMO {SHIFTS0[k]:.4f} -> {sh[k]:.4f}  "
              f"({d_ppm*reference_mhz:+.2f} Hz){flag}")
    print("\n===== coupling deviations from GISSMO (Hz) =====")
    cp = sol.x[N_SPINS:N_SPINS + len(CIDX)]
    for (i, j), c in zip(CIDX, cp):
        d = c - J0[i, j]
        flag = "  <-- moved" if abs(d) > float(CONFIG["reporting"]["coupling_flag_hz"]) else ""
        print(f"  J[{i:2d},{j:2d}] GISSMO {J0[i,j]:+.3f} -> {c:+.3f}  ({d:+.2f} Hz){flag}")

    # save the refined matrix
    out = os.path.join(HERE, "sucrose_matrix_MAP_refined.txt")
    np.savetxt(out, M_fit, fmt="%.9f", delimiter="\t")
    print(f"\nWrote refined matrix -> {out}")
    print("NEXT: the 'moved' params are only CANDIDATES. Check each one's")
    print("deviation reproduces across fields before believing it -- same logic")
    print("that rejected the CH2 splitting. A big deviation that only helps one")
    print("field is crowded-region overfitting, not a real refinement.")


if __name__ == "__main__":
    main()
