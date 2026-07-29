#!/usr/bin/env python3
"""Format the latest Bruker-style noesypr1d alanine fits.

The noesypr1d MATLAB drivers supply the experimental and fitted Spinach
traces.  The analytical AX3 curve is evaluated here from the same fixed
shifts/J matrix and the fitted field-specific linewidth.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs" / "alanine" / "georgia_noesypr1d_plots"
OUT.mkdir(parents=True, exist_ok=True)

HEDGES = "#B4BD00"       # experiment
GLORY_GLORY = "#E4002B"  # Spinach noesypr1d fit
OLYMPIC = "#004E60"      # analytical AX3 theory
BLACK = "#000000"
WHITE = "#FFFFFF"
GRID = "#D9E1E3"

SHIFT_HA = 3.7680
SHIFT_CH3 = 1.4655
J_HZ = 7.234

FIELDS = {
    600: (ROOT / "outputs" / "alanine" / "alanine_noesypr1d_experiment_fit_overlay.csv",
          ROOT / "outputs" / "alanine" / "alanine_noesypr1d_experiment_fit_summary.csv",
          599.764818881),
    900: (ROOT / "outputs" / "alanine" / "900MHz" / "alanine_noesypr1d_experiment_fit_overlay.csv",
          ROOT / "outputs" / "alanine" / "900MHz" / "alanine_noesypr1d_experiment_fit_summary.csv",
          899.794229013),
}


def read_csv(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ppm = np.array([float(row["ppm"]) for row in rows])
    exp = np.array([float(row["experiment_norm_interp"]) for row in rows])
    fit = np.array([float(row["noesypr1d_fit"]) for row in rows])
    return ppm, exp, fit


def read_summary(path: Path):
    with path.open(newline="") as handle:
        row = next(csv.DictReader(handle))
    return {key: float(value) for key, value in row.items()}


def theory_curve(ppm, sfo_mhz, lb_hz):
    j_ppm = J_HZ / sfo_mhz
    hwhm = (lb_hz / sfo_mhz) / 2.0

    def lorentz(x, center):
        return hwhm**2 / ((x - center) ** 2 + hwhm**2)

    curve = 1.5 * lorentz(ppm, SHIFT_CH3 - j_ppm / 2)
    curve += 1.5 * lorentz(ppm, SHIFT_CH3 + j_ppm / 2)
    centers = [SHIFT_HA - 1.5 * j_ppm, SHIFT_HA - 0.5 * j_ppm,
               SHIFT_HA + 0.5 * j_ppm, SHIFT_HA + 1.5 * j_ppm]
    weights = [1 / 8, 3 / 8, 3 / 8, 1 / 8]
    for center, weight in zip(centers, weights):
        curve += weight * lorentz(ppm, center)
    ch3 = (ppm >= 1.30) & (ppm <= 1.65)
    scale = np.max(curve[ch3])
    return curve / scale if scale > 0 else curve


def corr(a, b, mask):
    good = mask & np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[good], b[good])[0, 1])


def style(ax):
    ax.set_facecolor(WHITE)
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.8)
    ax.tick_params(axis="both", colors=BLACK, labelsize=15, width=1.4,
                   length=6)
    for spine in ax.spines.values():
        spine.set_color(BLACK)
        spine.set_linewidth(1.5)


def lines(ax, ppm, exp, fit, theory):
    h1, = ax.plot(ppm, exp, color=HEDGES, lw=2.5, label="Experiment")
    h2, = ax.plot(ppm, fit, color=GLORY_GLORY, lw=2.3, ls="--",
                  label="Spinach noesypr1d fit")
    h3, = ax.plot(ppm, theory, color=OLYMPIC, lw=2.2, ls=":",
                  label="AX₃ analytical theory")
    return [h1, h2, h3]


def make_field(field):
    overlay_path, summary_path, sfo = FIELDS[field]
    ppm, exp, fit = read_csv(overlay_path)
    summary = read_summary(summary_path)
    lb = summary["lb_Hz"]
    theory = theory_curve(ppm, sfo, lb)
    fit_mask = (((ppm >= 1.30) & (ppm <= 1.65)) |
                ((ppm >= 3.55) & (ppm <= 3.95)))
    r_fit = summary["r_noesypr1d_vs_expt"]
    r_theory = corr(theory, exp, fit_mask)

    full = (ppm >= 0.80) & (ppm <= 4.50)
    ymax = max(1.08, np.max(np.concatenate((exp[full], fit[full], theory[full]))) * 1.10)
    fig, ax = plt.subplots(figsize=(15.5, 8.5), dpi=180)
    fig.patch.set_facecolor(WHITE)
    handles = lines(ax, ppm[full], exp[full], fit[full], theory[full])
    style(ax)
    ax.set_xlim(4.50, 0.80)
    ax.set_ylim(-0.08, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=20,
                  fontweight="bold", color=BLACK)
    ax.set_ylabel("Normalised intensity", fontsize=20,
                  fontweight="bold", color=BLACK)
    ax.set_title(
        f"Alanine {field} MHz — Bruker-style noesypr1d validation\n"
        f"Experiment | Spinach noesypr1d fit | AX₃ theory   ·   "
        f"J = {J_HZ:.3f} Hz | LB = {lb:.3f} Hz | "
        f"r(Spinach, exp) = {r_fit:.4f} | r(AX₃, exp) = {r_theory:.4f}",
        fontsize=21, fontweight="bold", color=BLACK, pad=18)
    ax.legend(handles=handles, loc="upper left", fontsize=15, frameon=True,
              facecolor=WHITE, edgecolor=BLACK, framealpha=1.0)
    ax.text(0.99, 0.03,
            f"δHα = {SHIFT_HA:.4f} ppm   δCH₃ = {SHIFT_CH3:.4f} ppm\n"
            f"ppm offset = {summary['ppm_offset']:+.5f} ppm   "
            f"receiver phase = {summary['phase_deg']:+.1f}°",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=13,
            color=BLACK, bbox=dict(facecolor=WHITE, edgecolor=BLACK,
                                    linewidth=0.8, alpha=0.95))
    fig.tight_layout()
    fig.savefig(OUT / f"alanine_{field}MHz_noesypr1d_georgia_full.png",
                dpi=300, facecolor=WHITE, bbox_inches="tight")
    fig.savefig(OUT / f"alanine_{field}MHz_noesypr1d_georgia_full.pdf",
                facecolor=WHITE, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.0), dpi=180,
                             constrained_layout=True)
    fig.patch.set_facecolor(WHITE)
    for ax, lo, hi, title, upper in [
        (axes[0], 1.30, 1.65, "CH₃ doublet", 1.20),
        (axes[1], 3.55, 3.95, "Hα quartet", 0.55),
    ]:
        m = (ppm >= lo) & (ppm <= hi)
        handles = lines(ax, ppm[m], exp[m], fit[m], theory[m])
        style(ax)
        ax.set_xlim(hi, lo)
        ax.set_ylim(-0.05, upper)
        ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=17,
                      fontweight="bold", color=BLACK)
        ax.set_ylabel("Normalised intensity", fontsize=17,
                      fontweight="bold", color=BLACK)
        ax.set_title(title, fontsize=20, fontweight="bold", color=BLACK)
        ax.legend(handles=handles, loc="upper left", fontsize=11,
                  frameon=True, facecolor=WHITE, edgecolor=BLACK,
                  framealpha=1.0)
    axes[0].text(0.98, 0.93, f"δCH₃ = {SHIFT_CH3:.4f} ppm\nJ = {J_HZ:.3f} Hz",
                 transform=axes[0].transAxes, ha="right", va="top", fontsize=12,
                 color=BLACK, bbox=dict(facecolor=WHITE, edgecolor=BLACK,
                                        linewidth=0.8, alpha=0.95))
    axes[1].text(0.98, 0.93, f"δHα = {SHIFT_HA:.4f} ppm\nJ = {J_HZ:.3f} Hz",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=12,
                 color=BLACK, bbox=dict(facecolor=WHITE, edgecolor=BLACK,
                                        linewidth=0.8, alpha=0.95))
    fig.suptitle(
        f"Alanine {field} MHz — noesypr1d AX₃ multiplets\n"
        f"Spinach vs experiment r = {r_fit:.4f}   |   AX₃ vs experiment r = {r_theory:.4f}",
        fontsize=22, fontweight="bold", color=BLACK)
    fig.savefig(OUT / f"alanine_{field}MHz_noesypr1d_georgia_zooms.png",
                dpi=300, facecolor=WHITE, bbox_inches="tight")
    fig.savefig(OUT / f"alanine_{field}MHz_noesypr1d_georgia_zooms.pdf",
                facecolor=WHITE, bbox_inches="tight")
    plt.close(fig)
    print(f"{field} MHz: Spinach r={r_fit:.4f}; AX3 r={r_theory:.4f}; LB={lb:.4f} Hz")


def main():
    for field in (600, 900):
        make_field(field)
    print(f"Wrote noesypr1d Georgia plots to {OUT}")


if __name__ == "__main__":
    main()
