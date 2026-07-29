#!/usr/bin/env python3
"""Make Georgia-style plots for Mystery Sugar 1 (D-xylose candidate).

The source fit is the saved alpha/beta Spinach noesypr1d candidate from the
original Mystery Sugar 1 analysis.  Xylose does not have a deposited GISSMO
matrix, so the third trace is explicitly labelled as the equal-response
combined candidate model rather than GISSMO.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
PROJECT = ROOT.parent / "mystery_sugar_1_fit"
PREPARED = PROJECT / "prepared"
OVERLAY = PROJECT / "alpha_h4_structural_fit/candidate_h4_targeted/spinach_three_field_overlay.tsv"
SUMMARY = PROJECT / "alpha_h4_structural_fit/candidate_h4_targeted/spinach_validation_summary.txt"
OUT = ROOT / "outputs" / "mystery_sugar" / "georgia_plots"
OUT.mkdir(parents=True, exist_ok=True)

FIELDS = ("600", "900", "1100")
WINDOWS = {
    "alpha_anomeric": (5.30, 5.05),
    "beta_anomeric": (4.75, 4.40),
    "carbohydrate_region": (4.15, 3.10),
}
FULL = (5.35, 3.05)

# University of Georgia palette.
GLORY_GLORY = "#E4002B"  # refined alpha/beta Spinach candidate
HEDGES = "#B4BD00"       # experimental spectrum
OLYMPIC = "#004E60"      # equal-response candidate baseline
ARCH_BLACK = "#000000"
CHAPEL_BELL = "#FFFFFF"
GRID = "#D9E1E3"


def read_metrics() -> dict[str, dict[str, float]]:
    text = SUMMARY.read_text()
    metrics: dict[str, dict[str, float]] = {}
    pat = re.compile(
        r"(?m)^(600|900|1100) MHz\s+"
        r"alpha_lw_hz\s+([0-9.]+)\s+beta_lw_hz\s+([0-9.]+)\s+"
        r"combined_lw_hz\s+([0-9.]+)\s+"
        r"mixture_r\s+([0-9.]+)\s+mixture_rmse\s+([0-9.]+)\s+"
        r"combined_r\s+([0-9.]+)\s+combined_rmse\s+([0-9.]+)"
    )
    for m in pat.finditer(text):
        metrics[m.group(1)] = {
            "alpha_lw": float(m.group(2)),
            "beta_lw": float(m.group(3)),
            "combined_lw": float(m.group(4)),
            "candidate_r": float(m.group(5)),
            "candidate_rmse": float(m.group(6)),
            "baseline_r": float(m.group(7)),
            "baseline_rmse": float(m.group(8)),
        }
    return metrics


def read_field(field: str, overlay: pd.DataFrame) -> dict[str, pd.Series | float]:
    exp = pd.read_csv(PREPARED / f"mystery_sugar_1_{field}MHz_dss_aligned.tsv",
                      sep="\t")
    ov = overlay[overlay["field"].eq(f"{field} MHz")].copy()
    return {
        "exp": exp,
        "ppm": ov["ppm"],
        "candidate": ov["spinach_refined_mixture"],
        "baseline": ov["spinach_combined_matrix"],
    }


def style_axis(ax, fontsize: int = 13) -> None:
    ax.set_facecolor(CHAPEL_BELL)
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=fontsize, colors=ARCH_BLACK,
                   width=1.3, length=6)
    for spine in ax.spines.values():
        spine.set_color(ARCH_BLACK)
        spine.set_linewidth(1.5)


def draw(ax, data: dict[str, pd.Series | float], left: float, right: float):
    exp = data["exp"]
    ppm = data["ppm"]
    m_exp = (exp["ppm"] >= right) & (exp["ppm"] <= left)
    m_ov = (ppm >= right) & (ppm <= left)
    h_exp, = ax.plot(exp.loc[m_exp, "ppm"], exp.loc[m_exp, "intensity_norm"],
                     color=HEDGES, linewidth=2.2, label="Experimental spectrum",
                     zorder=3)
    h_cand, = ax.plot(ppm[m_ov], data["candidate"][m_ov], color=GLORY_GLORY,
                      linewidth=2.0, linestyle="--",
                      label="Refined α/β Spinach candidate", zorder=4)
    h_base, = ax.plot(ppm[m_ov], data["baseline"][m_ov], color=OLYMPIC,
                      linewidth=1.8, linestyle=":",
                      label="Equal-response candidate baseline", zorder=2)
    ymax = max(float(exp.loc[m_exp, "intensity_norm"].max()),
               float(data["candidate"][m_ov].max()),
               float(data["baseline"][m_ov].max()))
    return (h_exp, h_cand, h_base), max(1.02, 1.10 * ymax)


def save_individuals(fields: dict[str, dict], metrics: dict[str, dict]) -> None:
    for field in FIELDS:
        data = fields[field]
        m = metrics[field]
        # Full spectrum.
        fig, ax = plt.subplots(figsize=(14.5, 7.7), dpi=180)
        fig.patch.set_facecolor(CHAPEL_BELL)
        handles, ymax = draw(ax, data, *FULL)
        style_axis(ax, 14)
        ax.set_xlim(FULL)
        ax.set_ylim(-0.04, ymax)
        ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=18,
                      fontweight="bold")
        ax.set_ylabel("Normalised intensity", fontsize=18, fontweight="bold")
        ax.set_title(
            f"Mystery Sugar 1 — {field} MHz noesypr1d\n"
            f"D-xylose candidate | candidate r = {m['candidate_r']:.4f}, "
            f"RMSE = {m['candidate_rmse']:.4f}",
            fontsize=21, fontweight="bold", pad=16,
        )
        ax.legend(handles=handles, loc="upper left", fontsize=12,
                  frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                  framealpha=1.0)
        ax.text(0.99, 0.03,
                "Xylose has no deposited GISSMO matrix\n"
                "candidate baseline = equal-response α/β model",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
                bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                          linewidth=0.8, alpha=0.95))
        fig.tight_layout()
        fig.savefig(OUT / f"mystery_sugar1_xylose_{field}MHz_georgia_full.png",
                    dpi=300, bbox_inches="tight", facecolor=CHAPEL_BELL)
        fig.savefig(OUT / f"mystery_sugar1_xylose_{field}MHz_georgia_full.pdf",
                    bbox_inches="tight", facecolor=CHAPEL_BELL)
        plt.close(fig)

        # Three diagnostic zooms.
        fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.2), dpi=180,
                                 constrained_layout=True)
        fig.patch.set_facecolor(CHAPEL_BELL)
        handles = None
        for ax, (name, (left, right)) in zip(axes, WINDOWS.items()):
            handles, ymax = draw(ax, data, left, right)
            style_axis(ax, 12)
            ax.set_xlim(left, right)
            ax.set_ylim(-0.04, ymax)
            ax.set_title(name.replace("_", " "), fontsize=15, fontweight="bold")
            ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=12,
                          fontweight="bold")
            ax.set_ylabel("Normalised intensity", fontsize=12,
                          fontweight="bold")
        fig.suptitle(
            f"Mystery Sugar 1 — {field} MHz D-xylose candidate multiplet zooms",
            fontsize=21, fontweight="bold", color=ARCH_BLACK,
        )
        fig.legend(handles=handles, loc="upper center", ncol=3,
                   bbox_to_anchor=(0.5, 0.94), fontsize=11,
                   frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                   framealpha=1.0)
        fig.savefig(OUT / f"mystery_sugar1_xylose_{field}MHz_georgia_zooms.png",
                    dpi=300, bbox_inches="tight", facecolor=CHAPEL_BELL)
        fig.savefig(OUT / f"mystery_sugar1_xylose_{field}MHz_georgia_zooms.pdf",
                    bbox_inches="tight", facecolor=CHAPEL_BELL)
        plt.close(fig)


def save_combined(fields: dict[str, dict], metrics: dict[str, dict]) -> None:
    # Three-panel full-spectrum poster figure.
    fig, axes = plt.subplots(3, 1, figsize=(16.5, 12.0), dpi=180,
                             sharex=True)
    fig.patch.set_facecolor(CHAPEL_BELL)
    handles = None
    for ax, field in zip(axes, FIELDS):
        handles, ymax = draw(ax, fields[field], *FULL)
        style_axis(ax, 13)
        ax.set_xlim(FULL)
        ax.set_ylim(-0.04, ymax)
        m = metrics[field]
        ax.set_ylabel(f"{field} MHz\nNormalised", fontsize=13,
                      fontweight="bold")
        ax.text(0.01, 0.86,
                f"{field} MHz | candidate r = {m['candidate_r']:.4f} | "
                f"baseline r = {m['baseline_r']:.4f}",
                transform=ax.transAxes, ha="left", va="top", fontsize=13,
                fontweight="bold",
                bbox=dict(facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
                          linewidth=0.8, alpha=0.95))
    axes[-1].set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=17,
                        fontweight="bold")
    fig.suptitle("Mystery Sugar 1 — D-xylose candidate, three-field validation",
                 fontsize=23, fontweight="bold", y=0.98)
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.945), fontsize=12,
               frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
               framealpha=1.0)
    fig.subplots_adjust(top=0.88, bottom=0.08, left=0.10, right=0.98,
                        hspace=0.12)
    fig.savefig(OUT / "mystery_sugar1_xylose_georgia_full_3panel.png",
                dpi=300, bbox_inches="tight", facecolor=CHAPEL_BELL)
    fig.savefig(OUT / "mystery_sugar1_xylose_georgia_full_3panel.pdf",
                bbox_inches="tight", facecolor=CHAPEL_BELL)
    plt.close(fig)

    # Three fields by three scientifically useful regions.
    fig, axes = plt.subplots(3, 3, figsize=(18.5, 13.0), dpi=180,
                             constrained_layout=False)
    fig.patch.set_facecolor(CHAPEL_BELL)
    handles = None
    for row, field in enumerate(FIELDS):
        for col, (name, (left, right)) in enumerate(WINDOWS.items()):
            handles, ymax = draw(axes[row, col], fields[field], left, right)
            style_axis(axes[row, col], 11)
            axes[row, col].set_xlim(left, right)
            axes[row, col].set_ylim(-0.04, ymax)
            axes[row, col].set_title(f"{field} MHz — {name.replace('_', ' ')}",
                                     fontsize=13, fontweight="bold")
            axes[row, col].set_xlabel(r"$^1$H shift (ppm)", fontsize=11,
                                      fontweight="bold")
            if col == 0:
                axes[row, col].set_ylabel("Normalised", fontsize=11,
                                          fontweight="bold")
    fig.suptitle("Mystery Sugar 1 — D-xylose candidate multiplet zooms",
                 fontsize=23, fontweight="bold", y=0.98)
    fig.legend(handles=handles, loc="upper center", ncol=3,
               bbox_to_anchor=(0.5, 0.945), fontsize=11,
               frameon=True, facecolor=CHAPEL_BELL, edgecolor=ARCH_BLACK,
               framealpha=1.0)
    fig.subplots_adjust(top=0.88, bottom=0.07, left=0.06, right=0.99,
                        hspace=0.30, wspace=0.14)
    fig.savefig(OUT / "mystery_sugar1_xylose_georgia_zooms_3x3.png",
                dpi=300, bbox_inches="tight", facecolor=CHAPEL_BELL)
    fig.savefig(OUT / "mystery_sugar1_xylose_georgia_zooms_3x3.pdf",
                bbox_inches="tight", facecolor=CHAPEL_BELL)
    plt.close(fig)


def main() -> None:
    metrics = read_metrics()
    overlay = pd.read_csv(OVERLAY, sep="\t")
    fields = {field: read_field(field, overlay) for field in FIELDS}
    save_individuals(fields, metrics)
    save_combined(fields, metrics)
    print(f"Wrote Georgia-style Mystery Sugar 1 plots to {OUT}")


if __name__ == "__main__":
    main()
