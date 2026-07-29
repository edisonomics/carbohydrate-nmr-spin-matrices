#!/usr/bin/env python3
"""Make presentation-ready alanine AX3 plots from the latest Spinach traces.

The MATLAB drivers write one CSV per field containing the experimental trace,
the direct Spinach FFT, and the analytical first-order AX3 theory.  This
script only formats those already-generated traces; it does not refit them.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
INPUT_DIR = ROOT / "outputs" / "alanine"
OUTPUT_DIR = INPUT_DIR / "georgia_plots"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# University of Georgia palette requested for the poster/talk.
GLORY_GLORY = "#E4002B"  # Spinach candidate
HEDGES = "#B4BD00"       # experimental spectrum
OLYMPIC = "#004E60"       # analytical AX3 theory
ARCH_BLACK = "#000000"
CHAPEL_BELL = "#FFFFFF"
GRID = "#D9E1E3"

SHIFT_HA = 3.7680
SHIFT_CH3 = 1.4655
J_HZ = 7.234


def read_overlay(field: int):
    path = INPUT_DIR / f"alanine_spinach_fft_overlay_{field}MHz.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing latest overlay: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    ppm = np.array([float(row["ppm"]) for row in rows])
    spinach = np.array([float(row["spinach_norm"]) for row in rows])
    theory = np.array([float(row["theory_norm"]) for row in rows])
    experiment = np.array([float(row["experiment_interp_norm"]) for row in rows])
    return ppm, experiment, spinach, theory


def mask_for_fit(ppm: np.ndarray) -> np.ndarray:
    return (((ppm >= 1.30) & (ppm <= 1.65)) |
            ((ppm >= 3.55) & (ppm <= 3.95)))


def corr_rmse(a: np.ndarray, b: np.ndarray, mask: np.ndarray):
    good = mask & np.isfinite(a) & np.isfinite(b)
    aa, bb = a[good], b[good]
    return float(np.corrcoef(aa, bb)[0, 1]), float(np.sqrt(np.mean((aa - bb) ** 2)))


def style_axis(ax, fontsize=15):
    ax.set_facecolor(CHAPEL_BELL)
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.8)
    ax.tick_params(axis="both", labelsize=fontsize, colors=ARCH_BLACK,
                   width=1.4, length=6)
    for spine in ax.spines.values():
        spine.set_color(ARCH_BLACK)
        spine.set_linewidth(1.5)


def add_lines(ax, ppm, experiment, spinach, theory, label_prefix=""):
    # Draw the experiment first so the model traces can be compared directly.
    h_exp, = ax.plot(ppm, experiment, color=HEDGES, linewidth=2.5,
                     label="Experimental spectrum")
    h_sp, = ax.plot(ppm, spinach, color=GLORY_GLORY, linewidth=2.2,
                    linestyle="--", label="Spinach FFT")
    h_th, = ax.plot(ppm, theory, color=OLYMPIC, linewidth=2.2,
                    linestyle=":", label="AX₃ analytical theory")
    return h_exp, h_sp, h_th


def save_field_plots(field: int):
    ppm, experiment, spinach, theory = read_overlay(field)
    fit_mask = mask_for_fit(ppm)
    r_sp, rmse_sp = corr_rmse(spinach, experiment, fit_mask)
    r_th, rmse_th = corr_rmse(theory, experiment, fit_mask)
    r_st, rmse_st = corr_rmse(spinach, theory, fit_mask)

    # The current MATLAB output has tiny interpolation tails outside the
    # alanine window.  The publication view is intentionally restricted to
    # the actual alanine region.
    full = (ppm >= 0.80) & (ppm <= 4.50)
    y_full = np.concatenate((experiment[full], spinach[full], theory[full]))
    y_full = y_full[np.isfinite(y_full)]
    y_top = max(1.08, float(np.max(y_full)) * 1.10)

    # Full-spectrum figure.
    fig, ax = plt.subplots(figsize=(15.5, 8.5), dpi=180)
    fig.patch.set_facecolor(CHAPEL_BELL)
    handles = add_lines(ax, ppm[full], experiment[full], spinach[full], theory[full])
    style_axis(ax, 16)
    ax.set_xlim(4.50, 0.80)
    ax.set_ylim(-0.08, y_top)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=20,
                  fontweight="bold", color=ARCH_BLACK)
    ax.set_ylabel("Normalised intensity", fontsize=20,
                  fontweight="bold", color=ARCH_BLACK)
    ax.set_title(
        f"Alanine {field} MHz — experiment, Spinach FFT, and AX₃ theory\n"
        f"J = {J_HZ:.3f} Hz  |  r(Spinach, exp) = {r_sp:.4f}  |  "
        f"r(AX₃, exp) = {r_th:.4f}  |  RMSE = {rmse_sp:.4f}",
        fontsize=22, fontweight="bold", color=ARCH_BLACK, pad=18,
    )
    ax.legend(handles=handles, loc="upper left", fontsize=15,
              frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
              framealpha=1.0)
    ax.text(0.99, 0.03,
            f"δHα = {SHIFT_HA:.4f} ppm   δCH₃ = {SHIFT_CH3:.4f} ppm\n"
            f"Spinach vs AX₃: r = {r_st:.4f}, RMSE = {rmse_st:.4f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=13,
            color=ARCH_BLACK,
            bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                      linewidth=0.8, alpha=0.95))
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f"alanine_{field}MHz_georgia_full.png",
                dpi=300, facecolor=CHAPEL_BELL, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"alanine_{field}MHz_georgia_full.pdf",
                facecolor=CHAPEL_BELL, bbox_inches="tight")
    plt.close(fig)

    # Two diagnostic zooms: the AX3 methyl doublet and Halpha quartet.
    fig, axes = plt.subplots(1, 2, figsize=(16.5, 7.0), dpi=180,
                             constrained_layout=True)
    fig.patch.set_facecolor(CHAPEL_BELL)
    windows = [(1.30, 1.65, "CH₃ doublet", 1.20),
               (3.55, 3.95, "Hα quartet", 0.55)]
    for ax, (lo, hi, title, ymax) in zip(axes, windows):
        m = (ppm >= lo) & (ppm <= hi)
        handles = add_lines(ax, ppm[m], experiment[m], spinach[m], theory[m])
        style_axis(ax, 15)
        ax.set_xlim(hi, lo)
        ax.set_ylim(-0.05, ymax)
        ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=17,
                      fontweight="bold", color=ARCH_BLACK)
        ax.set_ylabel("Normalised intensity", fontsize=17,
                      fontweight="bold", color=ARCH_BLACK)
        ax.set_title(title, fontsize=20, fontweight="bold", color=ARCH_BLACK)
        ax.legend(handles=handles, loc="upper left", fontsize=11,
                  frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                  framealpha=1.0)
    axes[0].text(0.98, 0.93,
                 f"δCH₃ = {SHIFT_CH3:.4f} ppm\nJ = {J_HZ:.3f} Hz",
                 transform=axes[0].transAxes, ha="right", va="top",
                 fontsize=12, color=ARCH_BLACK,
                 bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                           linewidth=0.8, alpha=0.95))
    axes[1].text(0.98, 0.93,
                 f"δα = {SHIFT_HA:.4f} ppm\nJ = {J_HZ:.3f} Hz",
                 transform=axes[1].transAxes, ha="right", va="top",
                 fontsize=12, color=ARCH_BLACK,
                 bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                           linewidth=0.8, alpha=0.95))
    fig.suptitle(
        f"Alanine {field} MHz — AX₃ multiplet detail\n"
        f"Spinach vs experiment r = {r_sp:.4f}  |  AX₃ vs experiment r = {r_th:.4f}",
        fontsize=22, fontweight="bold", color=ARCH_BLACK,
    )
    fig.savefig(OUTPUT_DIR / f"alanine_{field}MHz_georgia_zooms.png",
                dpi=300, facecolor=CHAPEL_BELL, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"alanine_{field}MHz_georgia_zooms.pdf",
                facecolor=CHAPEL_BELL, bbox_inches="tight")
    plt.close(fig)

    return dict(field=field, r_spinach=r_sp, rmse_spinach=rmse_sp,
                r_theory=r_th, rmse_theory=rmse_th,
                r_spinach_theory=r_st, rmse_spinach_theory=rmse_st)


def save_combined_plots():
    """Write poster-friendly 600/900 MHz comparison figures.

    The first figure is a two-panel full-spectrum comparison.  The second is
    a four-panel zoom figure (field strength by multiplet), so every panel has
    enough horizontal space for the three traces and its legend.
    """
    traces = {}
    for field in (600, 900):
        ppm, experiment, spinach, theory = read_overlay(field)
        fit_mask = mask_for_fit(ppm)
        r_sp, rmse_sp = corr_rmse(spinach, experiment, fit_mask)
        r_th, rmse_th = corr_rmse(theory, experiment, fit_mask)
        traces[field] = dict(ppm=ppm, experiment=experiment,
                             spinach=spinach, theory=theory,
                             r_sp=r_sp, rmse_sp=rmse_sp,
                             r_th=r_th, rmse_th=rmse_th)

    # Two-panel full spectra, with one shared legend.
    fig, axes = plt.subplots(2, 1, figsize=(16.5, 10.5), dpi=180,
                             sharex=True, constrained_layout=False)
    fig.patch.set_facecolor(CHAPEL_BELL)
    handles = None
    for ax, field in zip(axes, (600, 900)):
        data = traces[field]
        full = (data["ppm"] >= 0.80) & (data["ppm"] <= 4.50)
        handles = add_lines(ax, data["ppm"][full], data["experiment"][full],
                            data["spinach"][full], data["theory"][full])
        style_axis(ax, 15)
        y = np.concatenate((data["experiment"][full], data["spinach"][full],
                            data["theory"][full]))
        ax.set_xlim(4.50, 0.80)
        ax.set_ylim(-0.06, max(1.08, float(np.nanmax(y)) * 1.10))
        ax.set_ylabel(f"{field} MHz\nNormalised intensity", fontsize=15,
                      fontweight="bold", color=ARCH_BLACK)
        ax.text(0.01, 0.88,
                f"{field} MHz  |  Spinach r = {data['r_sp']:.4f}  |  "
                f"AX₃ r = {data['r_th']:.4f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=15,
                fontweight="bold", color=ARCH_BLACK,
                bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                          linewidth=0.8, alpha=0.95))
    axes[-1].set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=19,
                        fontweight="bold", color=ARCH_BLACK)
    fig.suptitle("Alanine — two-field full-spectrum validation",
                 fontsize=24, fontweight="bold", color=ARCH_BLACK, y=0.98)
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.94), fontsize=14,
               frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
               framealpha=1.0)
    fig.subplots_adjust(top=0.86, bottom=0.10, left=0.10, right=0.98,
                        hspace=0.10)
    fig.savefig(OUTPUT_DIR / "alanine_600_900MHz_georgia_full_2panel.png",
                dpi=300, facecolor=CHAPEL_BELL, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "alanine_600_900MHz_georgia_full_2panel.pdf",
                facecolor=CHAPEL_BELL, bbox_inches="tight")
    plt.close(fig)

    # Four-panel zooms: rows are fields, columns are the two AX3 multiplets.
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 10.0), dpi=180,
                             constrained_layout=False)
    fig.patch.set_facecolor(CHAPEL_BELL)
    handles = None
    windows = [(1.30, 1.65, "CH₃ doublet", 1.20),
               (3.55, 3.95, "Hα quartet", 0.55)]
    for row, field in enumerate((600, 900)):
        data = traces[field]
        for col, (lo, hi, title, ymax) in enumerate(windows):
            ax = axes[row, col]
            m = (data["ppm"] >= lo) & (data["ppm"] <= hi)
            handles = add_lines(ax, data["ppm"][m], data["experiment"][m],
                                data["spinach"][m], data["theory"][m])
            style_axis(ax, 14)
            ax.set_xlim(hi, lo)
            ax.set_ylim(-0.04, ymax)
            ax.set_title(f"{field} MHz — {title}", fontsize=17,
                         fontweight="bold", color=ARCH_BLACK)
            ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=15,
                          fontweight="bold", color=ARCH_BLACK)
            if col == 0:
                ax.set_ylabel("Normalised intensity", fontsize=15,
                              fontweight="bold", color=ARCH_BLACK)
            ax.text(0.98, 0.93,
                    (f"δCH₃ = {SHIFT_CH3:.4f} ppm\nJ = {J_HZ:.3f} Hz"
                     if col == 0 else
                     f"δα = {SHIFT_HA:.4f} ppm\nJ = {J_HZ:.3f} Hz"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=11,
                    color=ARCH_BLACK,
                    bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                              linewidth=0.7, alpha=0.95))
    fig.suptitle("Alanine — AX₃ multiplet zooms at two fields",
                 fontsize=24, fontweight="bold", color=ARCH_BLACK, y=0.98)
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.94), fontsize=13,
               frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
               framealpha=1.0)
    fig.subplots_adjust(top=0.85, bottom=0.09, left=0.08, right=0.98,
                        hspace=0.28, wspace=0.15)
    fig.savefig(OUTPUT_DIR / "alanine_600_900MHz_georgia_zooms_2panel.png",
                dpi=300, facecolor=CHAPEL_BELL, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / "alanine_600_900MHz_georgia_zooms_2panel.pdf",
                facecolor=CHAPEL_BELL, bbox_inches="tight")
    plt.close(fig)


def main():
    results = [save_field_plots(field) for field in (600, 900)]
    save_combined_plots()
    print(f"Wrote Georgia-style alanine plots to {OUTPUT_DIR}")
    for result in results:
        print(result)


if __name__ == "__main__":
    main()
