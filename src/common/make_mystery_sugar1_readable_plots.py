#!/usr/bin/env python3
"""Make publication-readable Mystery Sugar 1 (xylose candidate) overlays.

The source data are the DSS-aligned experimental spectra and the numerical
Spinach overlay exported by the alpha/beta candidate fit.  This deliberately
does not relabel the equal-response model as GISSMO: xylose did not have a
deposited GISSMO matrix.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


REPO = Path("/Users/cece/Desktop/edison_lab/final_repo")
PROJECT = REPO.parent / "mystery_sugar_1_fit"
PREPARED = PROJECT / "prepared"
OVERLAY = PROJECT / "alpha_h4_structural_fit/candidate_h4_targeted/spinach_three_field_overlay.tsv"
SUMMARY = PROJECT / "alpha_h4_structural_fit/candidate_h4_targeted/spinach_validation_summary.txt"
OUT = PROJECT / "alpha_h4_structural_fit/candidate_h4_targeted/readable_plots"

FIELDS = ("600", "900", "1100")
WINDOWS = {
    "full": (5.35, 3.05),
    "alpha_anomeric": (5.30, 5.05),
    "beta_anomeric": (4.75, 4.40),
    "carbohydrate_region": (4.15, 3.10),
}


def read_metrics() -> dict[str, tuple[float, float]]:
    text = SUMMARY.read_text()
    metrics: dict[str, tuple[float, float]] = {}
    for m in re.finditer(
        r"(?m)^\s*(600|900|1100)\s+MHz.*?mixture_r\s+([0-9.]+)\s+mixture_rmse\s+([0-9.]+)",
        text,
    ):
        metrics[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return metrics


def style_axis(ax, left: float, right: float, title: str, ymax: float) -> None:
    ax.set_xlim(left, right)
    ax.set_ylim(-0.04 * ymax, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=13)
    ax.set_ylabel("Normalized intensity", fontsize=13)
    ax.set_title(title, fontsize=15, weight="bold", pad=10)
    ax.grid(True, color="#d9d9d9", linewidth=0.8, alpha=0.65)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        spine.set_color("#222222")
    ax.tick_params(labelsize=11, width=1.2)


def plot_one(field: str, window_name: str, left: float, right: float, metrics) -> Path:
    exp_path = PREPARED / f"mystery_sugar_1_{field}MHz_dss_aligned.tsv"
    exp = pd.read_csv(exp_path, sep="\t")
    ov = pd.read_csv(OVERLAY, sep="\t")
    ov = ov[ov["field"].eq(f"{field} MHz")]

    exp = exp[(exp.ppm >= right) & (exp.ppm <= left)]
    ov = ov[(ov.ppm >= right) & (ov.ppm <= left)]
    ymax = max(
        float(exp.intensity_norm.max()),
        float(ov.spinach_refined_mixture.max()),
        float(ov.spinach_combined_matrix.max()),
    ) * 1.12

    r, rmse = metrics.get(field, (float("nan"), float("nan")))
    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(exp.ppm, exp.intensity_norm, color="black", linewidth=1.8,
            label="Experimental spectrum", zorder=3)
    ax.plot(ov.ppm, ov.spinach_refined_mixture, color="#d62728", linewidth=2.2,
            linestyle="--", label="Refined α/β Spinach mixture", zorder=4)
    ax.plot(ov.ppm, ov.spinach_combined_matrix, color="#1769aa", linewidth=1.8,
            linestyle=":", label="Equal-response combined model", zorder=2)

    if window_name == "full":
        title = (f"Mystery Sugar 1 — {field} MHz | refined α/β mixture "
                 f"r = {r:.4f}, RMSE = {rmse:.4f}")
    else:
        title = f"Mystery Sugar 1 — {field} MHz | {window_name.replace('_', ' ')}"
    style_axis(ax, left, right, title, ymax)
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="#333333",
              fontsize=11)
    ax.text(
        0.99, 0.97,
        "black: experiment\nred dashed: refined α/β mixture\nblue dotted: equal-response model",
        transform=ax.transAxes, ha="right", va="top", fontsize=10,
        color="#222222", bbox=dict(facecolor="white", edgecolor="#888888", alpha=0.94),
    )
    fig.tight_layout()
    out = OUT / f"mystery_sugar_1_{field}MHz_{window_name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = read_metrics()
    made = []
    for field in FIELDS:
        for name, (left, right) in WINDOWS.items():
            made.append(plot_one(field, name, left, right, metrics))
    print(f"Wrote {len(made)} high-contrast plots to {OUT}")
    for p in made:
        print(p)


if __name__ == "__main__":
    main()
