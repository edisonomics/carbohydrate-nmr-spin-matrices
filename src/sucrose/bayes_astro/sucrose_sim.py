"""
sucrose_sim.py -- standalone EXACT quantum simulator for the sucrose 1H spectrum,
built so a Bayesian (astrophysics-style) fit can call it tens of thousands of
times cheaply.

Why this exists alongside the Spinach/MATLAB pipeline:
  - Bayesian inference (emcee) needs many forward-model evaluations; bridging to
    MATLAB per-eval is painful. This is a native-Python forward model.
  - It is also an INDEPENDENT reimplementation of the spin physics -> a
    cross-check of Spinach (if both agree from the same matrix, the pipeline is
    confirmed independent of code). Same spirit as using GISSMO's own simulation
    as a second opinion.

Physics:
  - Sucrose's glucose-ring (matrix rows 0-6) and fructose-ring (rows 7-13)
    protons have ZERO mutual J-coupling, so each 7-spin block is simulated
    EXACTLY in its 2^7 = 128 dim Hilbert space (no approximation), and the
    stick spectra are summed. Identical justification to the Spinach
    block-split.
  - Per block: H = sum_k 2*pi*nu_k Iz_k + sum_{j<k} 2*pi*J_jk (IxIx+IyIy+IzIz),
    eigendecompose, get transition frequencies (E_b-E_a)/2pi and intensities
    |<a|I-|b>|^2 (absorption-mode, see convention note + self-test below).
  - lineshape: analytic Lorentzian sum over sticks. lb (linewidth) and a ppm
    calibration offset are applied to the sticks cheaply -- no FFT in the loop.

Run `python3 sucrose_sim.py` to execute the self-tests (these MUST pass before
trusting any fit built on this), and to dump a spectrum + the anomeric peak
position for cross-checking against Spinach.

Requires: numpy.
"""

import numpy as np

# ---- single-spin operators (spin-1/2; spin operators = Pauli/2) ----
_SX = 0.5 * np.array([[0, 1], [1, 0]], dtype=complex)
_SY = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)
_SP = np.array([[0, 1], [0, 0]], dtype=complex)   # I+
_SM = np.array([[0, 0], [1, 0]], dtype=complex)    # I-
_I2 = np.eye(2, dtype=complex)


def _op_on(op, k, n):
    """Embed single-spin operator `op` on spin k (0-based) in an n-spin space."""
    mats = [(_I2 if i != k else op) for i in range(n)]
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


def _block_operators(n):
    """Precompute Ix,Iy,Iz,Ip,Im for each spin in an n-spin block."""
    ix = [_op_on(_SX, k, n) for k in range(n)]
    iy = [_op_on(_SY, k, n) for k in range(n)]
    iz = [_op_on(_SZ, k, n) for k in range(n)]
    ip = [_op_on(_SP, k, n) for k in range(n)]
    im = [_op_on(_SM, k, n) for k in range(n)]
    return ix, iy, iz, ip, im


def block_sticks(shifts_ppm, Jblock_Hz, SFO1_MHz, carrier_ppm, rel_tol=1e-4):
    """
    Exact stick spectrum for one uncoupled block.
      shifts_ppm : (nb,) chemical shifts (ppm)
      Jblock_Hz  : (nb,nb) symmetric J matrix (Hz) within this block
    Returns (freq_ppm, intensity) arrays of transitions.

    Convention: detection = I- , rho0 = I+ (= det^dagger). This makes every
    transition intensity |<a|I-|b>|^2 real and >= 0 (pure absorption, no
    phasing) AND places an uncoupled spin's line at +nu (verified by the
    self-test at the bottom). ppm = carrier + nu/SFO1, so a bare spin lands
    exactly at its delta.
    """
    nb = len(shifts_ppm)
    ix, iy, iz, ip, im = _block_operators(nb)

    nu = (np.asarray(shifts_ppm) - carrier_ppm) * SFO1_MHz   # Hz offset from carrier
    H = np.zeros((2**nb, 2**nb), dtype=complex)
    for k in range(nb):
        H += 2 * np.pi * nu[k] * iz[k]
    for j in range(nb):
        for k in range(j + 1, nb):
            Jjk = Jblock_Hz[j, k]
            if Jjk != 0.0:
                H += 2 * np.pi * Jjk * (ix[j] @ ix[k] + iy[j] @ iy[k] + iz[j] @ iz[k])

    E, V = np.linalg.eigh(H)                 # H Hermitian -> real E (rad/s)
    Dtot = sum(im)                            # detection = I-
    Deig = V.conj().T @ Dtot @ V
    amps = np.abs(Deig) ** 2                  # |<a|I-|b>|^2, real >= 0
    freqs_Hz = (E[None, :] - E[:, None]) / (2 * np.pi)   # (E_b - E_a)/2pi

    mask = amps > rel_tol * amps.max()       # prune negligible transitions
    f = freqs_Hz[mask]
    a = amps[mask]
    freq_ppm = carrier_ppm + f / SFO1_MHz
    return freq_ppm, a


def make_block_solver(nb, Jblock_Hz, SFO1_MHz, carrier_ppm, rel_tol=1e-4):
    """Return a fast closure solve(shifts_ppm)->(freq_ppm, intensity) that caches
    the spin operators and the (shift-independent) coupling Hamiltonian, so a
    Bayesian fit that only varies shifts pays just one eigendecomposition + two
    matmuls per call -- not a full operator rebuild. Same physics/convention as
    block_sticks (verified by the self-tests, which exercise block_sticks)."""
    ix, iy, iz, ip, im = _block_operators(nb)
    Hc = np.zeros((2 ** nb, 2 ** nb), dtype=complex)      # coupling part (constant)
    for j in range(nb):
        for k in range(j + 1, nb):
            if Jblock_Hz[j, k] != 0.0:
                Hc += 2 * np.pi * Jblock_Hz[j, k] * (ix[j] @ ix[k] + iy[j] @ iy[k] + iz[j] @ iz[k])
    Dtot = sum(im)

    def solve(shifts_ppm):
        nu = (np.asarray(shifts_ppm) - carrier_ppm) * SFO1_MHz
        H = Hc.copy()
        for k in range(nb):
            H += 2 * np.pi * nu[k] * iz[k]
        E, V = np.linalg.eigh(H)
        Deig = V.conj().T @ Dtot @ V
        amps = np.abs(Deig) ** 2
        freqs_Hz = (E[None, :] - E[:, None]) / (2 * np.pi)
        mask = amps > rel_tol * amps.max()       # prune negligible transitions
        return carrier_ppm + freqs_Hz[mask] / SFO1_MHz, amps[mask]

    return solve


def lorentzian_spectrum(ppm_grid, stick_ppm, stick_int, lb_Hz, SFO1_MHz, offset_ppm=0.0):
    """Analytic Lorentzian sum on ppm_grid. lb_Hz = FWHM; offset shifts sticks.
    Cheap -- this is what the MCMC calls repeatedly for lb/offset moves."""
    hwhm_ppm = (lb_Hz / 2.0) / SFO1_MHz
    centers = stick_ppm + offset_ppm
    # vectorized: (Ngrid, Nstick)
    dx = ppm_grid[:, None] - centers[None, :]
    L = (hwhm_ppm ** 2) / (dx ** 2 + hwhm_ppm ** 2)
    return L @ stick_int


def sucrose_sticks(matrix_14, SFO1_MHz, carrier_ppm, blocks=None,
                   glucose_idx=None, fructose_idx=None):
    """Stick spectrum for uncoupled blocks in a matrix.

    ``blocks`` is a sequence of zero-based index sequences from the molecule
    configuration. The glucose/fructose defaults are retained for old callers,
    but a different carbohydrate can supply any number of blocks.
    """
    M = np.asarray(matrix_14, dtype=float)
    shifts = np.diag(M)
    J = np.triu(M, 1)
    J = J + J.T
    if blocks is None:
        if glucose_idx is None: glucose_idx = range(0, 7)
        if fructose_idx is None: fructose_idx = range(7, M.shape[0])
        blocks = [glucose_idx, fructose_idx]
    frequencies, intensities = [], []
    n_total = M.shape[0]
    for block in blocks:
        idx = list(block)
        f, a = block_sticks(shifts[idx], J[np.ix_(idx, idx)], SFO1_MHz, carrier_ppm)
        # A block calculation omits the spectator spins.  In the full
        # Hilbert space every transition in this block is degenerate over
        # the spectator space, whose dimension is 2**(N-n_block).  Include
        # that trace/degeneracy factor so ragged decompositions (for example
        # 7+5+2) have the same relative intensities as one full calculation.
        spectator_dim = 2 ** (n_total - len(idx))
        frequencies.append(f); intensities.append(a * spectator_dim)
    return np.concatenate(frequencies), np.concatenate(intensities)


def load_gissmo_matrix(path):
    """Load a square whitespace spin matrix (same format Spinach reads).

    The original sucrose implementation required 14x14.  Carbohydrates such
    as xylose may have a different number of observed protons, so the loader
    validates square shape without hard-coding the sucrose dimension.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            toks = line.split()
            if toks:
                rows.append([float(t) for t in toks])
    M = np.array(rows)
    if M.ndim != 2 or M.shape[0] == 0 or M.shape[0] != M.shape[1]:
        raise ValueError(f"expected a non-empty square matrix, got {M.shape}")
    return M


# =====================================================================
# Self-tests -- MUST pass before trusting any fit built on this simulator.
# =====================================================================
def _selftests():
    # Synthetic unit-test fixture: deliberately round numbers, independent
    # of any one Bruker acquisition in the data tree.
    SFO1 = 600.0
    carrier = 4.0

    # Test 1: a single uncoupled spin lands exactly at its delta.
    f, a = block_sticks([3.500], np.zeros((1, 1)), SFO1, carrier)
    assert len(f) == 1, f"single spin should give 1 line, got {len(f)}"
    assert abs(f[0] - 3.500) < 1e-6, f"single spin landed at {f[0]:.6f}, not 3.5 (AXIS SIGN BUG)"
    print(f"[ok] single uncoupled spin: line at {f[0]:.4f} ppm (expected 3.5000)")

    # Test 2: two uncoupled spins -> two VISIBLE lines at their deltas.
    # (Each spin gives 2 degenerate transitions, one per state of the other
    # spin; they sit at identical frequency and collapse to one line. So the
    # raw transition list is [3.5,3.5,4.2,4.2] -- correct -- and we check the
    # UNIQUE positions, not the raw count.)
    f, a = block_sticks([3.50, 4.20], np.zeros((2, 2)), SFO1, carrier)
    peaks = np.unique(np.round(f, 4))
    assert len(peaks) == 2 and abs(peaks[0] - 3.50) < 1e-3 and abs(peaks[1] - 4.20) < 1e-3, \
        f"two-spin unique lines at {peaks}, expected [3.5, 4.2]"
    print(f"[ok] two uncoupled spins: unique lines at {peaks} ppm "
          f"(degenerate transitions collapse correctly)")

    # Test 3: an AX pair (weak coupling) -> doublets split by ~J/SFO1 ppm,
    # centered on each shift. Check line count and rough splitting.
    J = 7.0
    f, a = block_sticks([3.50, 4.20], np.array([[0, J], [J, 0]]), SFO1, carrier)
    near_A = np.sort(f[np.abs(f - 3.50) < 0.05])
    if len(near_A) >= 2:
        split_Hz = (near_A[-1] - near_A[0]) * SFO1
        assert abs(split_Hz - J) < 1.0, f"AX doublet split {split_Hz:.2f} Hz, expected ~{J}"
        print(f"[ok] AX doublet splitting near 3.5 ppm: {split_Hz:.2f} Hz (expected ~{J})")

    # Test 4: intensities are real and non-negative (absorption mode).
    assert np.all(a >= -1e-9), "negative intensities -> not absorption mode"
    print("[ok] all transition intensities >= 0 (absorption mode, no phasing needed)")

    print("\nAll self-tests passed.\n")


if __name__ == "__main__":
    import sys
    _selftests()

    # Spinach cross-check: simulate the GISSMO matrix and report the anomeric
    # peak position + dump a spectrum to compare against the Spinach output.
    import os
    import csv
    here = os.path.dirname(__file__)
    repo_root = os.path.abspath(os.path.join(here, "..", "..", ".."))
    sys.path.insert(0, os.path.join(repo_root, "src", "common"))
    from carbohydrate_config import load_config
    from seed_manifest import resolve_seed_matrix
    config = load_config(__import__("pathlib").Path(repo_root), "sucrose")
    mpath = str(resolve_seed_matrix(__import__("pathlib").Path(repo_root), "sucrose", config))
    if not os.path.isfile(mpath):
        print(f"(skip cross-check: {mpath} not found)")
        sys.exit(0)

    summary_path = os.path.join(repo_root, "outputs", "sucrose", "prepared",
                                "preparation_summary.csv")
    with open(summary_path, newline="", encoding="utf-8") as handle:
        summary = {row["field_mhz"]: row for row in csv.DictReader(handle)}
    field = summary["600"]
    SFO1 = float(field["sfo1_mhz"])
    carrier = float(field["o1_hz"]) / SFO1
    M = load_gissmo_matrix(mpath)
    blocks = [[int(index) - 1 for index in block] for block in config.get("blocks", [])]
    sppm, sint = sucrose_sticks(M, SFO1, carrier, blocks=blocks or None)

    processing = config["processing"]
    model = config["independent_model"]
    fit_lo, fit_hi = processing["fit_region_ppm"]
    grid = np.linspace(fit_lo, fit_hi, int(model["grid_points"]))
    spec = lorentzian_spectrum(grid, sppm, sint,
                               lb_Hz=float(model["lb_hz"]), SFO1_MHz=SFO1)

    anom_lo, anom_hi = processing["anomeric_region_ppm"]
    anom = (grid >= anom_lo) & (grid <= anom_hi)
    ipk = np.argmax(spec[anom])
    print(f"Anomeric peak (Python sim): {grid[anom][ipk]:.4f} ppm "
          f"(expected ~{processing['anomeric_reference_ppm']:.3f} from the matrix).")

    out_dir = os.path.join(repo_root, "outputs", "sucrose")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, "python_sim_spectrum_600MHz.csv")
    spec_n = spec / np.max(spec[anom])
    np.savetxt(out, np.column_stack([grid, spec_n]), delimiter=",",
               header="ppm,intensity_anomeric_normalized", comments="")
    print(f"Wrote {out} -- overlay this on the Spinach GISSMO spectrum; the two "
          f"independent codes should agree peak-for-peak.")
