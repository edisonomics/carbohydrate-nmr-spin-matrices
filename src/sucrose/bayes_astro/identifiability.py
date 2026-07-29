"""Identifiability analysis of the sucrose spin matrix under 1D data.

Builds the Jacobian J = d(spectrum)/d(parameters) at the GISSMO matrix, in
per-Hz units (so shifts and couplings are comparable), for single-field vs
multi-field data. SVD(J) then quantifies EXACTLY how non-unique the matrix is:

  - singular-value spectrum + condition number     (degeneracy severity)
  - effective rank vs the noise floor              (# of IDENTIFIABLE directions)
  - (n_params - rank) = # of degenerate directions (the non-uniqueness dimension)
  - CRLB per-parameter uncertainty (Hz)            (which params are undetermined)
  - the near-null right-singular vectors           (WHICH combos are unidentifiable)

Local identifiability <=> J full column rank. This is the mathematical test of
"is the matrix uniquely determined by this data."
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import joint_fit as jf

SHIFTS0 = jf.SHIFTS0
CIDX = jf.CIDX
coup0 = np.array([jf.J0[i, j] for (i, j) in CIDX])
NAMES = jf.NAMES
DATA = jf.DATA
NP = jf.N_SPINS + len(CIDX)
LABELS = [f"d{k}" for k in range(jf.N_SPINS)] + [f"J{i}{j}" for (i, j) in CIDX]
SFO1_REF = float(jf.CONFIG["identifiability"]["reference_field_mhz"])
initial_lb = list(jf.MODEL_CONFIG["initial_linewidth_hz"])
LB = dict(zip(NAMES, initial_lb))
SIG = float(jf.CONFIG["identifiability"]["noise_sigma"])


def _M(shifts, coups):
    M = np.diag(shifts).astype(float)
    for (i, j), c in zip(CIDX, coups):
        M[i, j] = c
    return M


def forward(shifts, coups, fields):
    M = _M(shifts, coups)
    return np.concatenate([jf.field_model(M, DATA[n], LB[n], 0.0) for n in fields])


def jac(fields):
    ds = float(jf.CONFIG["identifiability"]["shift_step_ppm"])
    dj = float(jf.CONFIG["identifiability"]["coupling_step_hz"])
    cols = []
    for k in range(jf.N_SPINS):
        sp = SHIFTS0.copy(); sp[k] += ds
        sm = SHIFTS0.copy(); sm[k] -= ds
        cols.append((forward(sp, coup0, fields) - forward(sm, coup0, fields)) / (2 * ds * SFO1_REF))
    for k in range(len(coup0)):
        cp = coup0.copy(); cp[k] += dj
        cm = coup0.copy(); cm[k] -= dj
        cols.append((forward(SHIFTS0, cp, fields) - forward(SHIFTS0, cm, fields)) / (2 * dj))
    return np.array(cols).T                    # (ndata, NP): d(spec)/d(param[Hz])


def analyze(tag, fields):
    J = jac(fields)
    ndata = J.shape[0]
    U, S, Vt = np.linalg.svd(J, full_matrices=False)
    noise_L2 = SIG * np.sqrt(ndata)            # a 1-Hz move is visible iff sv_i > noise_L2
    print(f"\n===== {tag}   (Jacobian {J.shape}, {NP} params) =====")
    print("singular values:\n  " + np.array2string(S, precision=2, max_line_width=150))
    print(f"condition number sv_max/sv_min = {S[0]/S[-1]:.2e}")
    print(f"noise floor (a 1-Hz change is detectable iff sv > {noise_L2:.3f}):")
    ident = int((S > noise_L2).sum())
    print(f"  IDENTIFIABLE directions (sv > noise floor): {ident} / {NP}"
          f"   -> {NP - ident} DEGENERATE direction(s)")
    for thr in (1e-2, 1e-3):
        print(f"  (rank at sv > {thr:g}*sv_max: {int((S > thr*S[0]).sum())}/{NP})")

    crlb = SIG * np.sqrt(Vt.T**2 @ (1.0 / S**2))
    order = np.argsort(-crlb)
    print("\n  LEAST determined parameters (CRLB, Hz; huge => unidentifiable):")
    for o in order[:8]:
        u = crlb[o]
        print(f"    {LABELS[o]:>5}: +/- {u:.1f} Hz" if u < 1e4 else f"    {LABELS[o]:>5}: +/- unbounded")
    print("  BEST determined:")
    for o in order[::-1][:4]:
        print(f"    {LABELS[o]:>5}: +/- {crlb[o]:.3f} Hz")

    print("\n  UNIDENTIFIABLE combinations (smallest singular values):")
    for idx in range(NP - 1, NP - 4, -1):
        v = Vt[idx]; top = np.argsort(-np.abs(v))[:5]
        print(f"    sv={S[idx]:.2e}:  " + "  ".join(f"{v[t]:+.2f}{LABELS[t]}" for t in top))
    return S


single_name = NAMES[0]
configured_name = " + ".join(NAMES)
S1 = analyze(f"SINGLE FIELD ({single_name})", [single_name])
S3 = analyze(f"CONFIGURED FIT FIELDS ({configured_name})", NAMES)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8.5, 5.2))
ax.semilogy(range(1, NP + 1), S1 / S1[0], "o-", label="600n only", color="tab:blue")
ax.semilogy(range(1, NP + 1), S3 / S3[0], "s-", label="all 3 fields", color="tab:red")
ax.axhline(1e-2, ls="--", color="0.5", lw=1.4, label="~noise floor")
ax.set_xlabel("singular-value index (1 = best-determined)", fontsize=12)
ax.set_ylabel("singular value / max", fontsize=12)
ax.set_title("Sucrose matrix identifiability under 1D data\n"
             "singular values below the floor = parameter directions the data can't determine",
             fontsize=12)
ax.legend(fontsize=11); ax.grid(alpha=0.3, which="both")
fig.tight_layout()
fig.savefig(os.path.join(os.path.dirname(__file__), "identifiability_svd.png"), dpi=160)
print("\nwrote identifiability_svd.png")
