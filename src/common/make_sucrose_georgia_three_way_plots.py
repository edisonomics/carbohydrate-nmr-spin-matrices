#!/usr/bin/env python3
"""Make readable sucrose experiment/Spinach/GISSMO figures.

The numerical curves are the archived three-way comparisons from the original
sucrose Spinach work.  This plotting layer changes only presentation: it does
not refit the matrix or alter any spectrum.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


ROOT = Path("/Users/cece/Desktop/edison_lab")
SOURCE = ROOT / "sucrose_spinach_comparison/spinach_curve_exports"
OUT = ROOT / "final_repo/outputs/sucrose/georgia_three_way_plots"

# Georgia palette: Hedges, Glory Glory, and Olympic. Arch Black remains the
# neutral for axes and text so all three data series are color-coded.
BLACK = "#000000"
HEDGES = "#B4BD00"
GLORY = "#E4002B"
OLYMPIC = "#004E60"
GRID = "#D9E2E5"

FIELDS = ("600", "900", "1100")
MODEL = "ch2_intrinsic_voigt_direct"

# Standard diagnostic windows and two narrower “super zooms” around the most
# informative anomeric/crowded features.
WINDOWS = {
    "full_sugar_region": (5.80, 3.00),
    "anomeric": (5.45, 5.35),
    "fructose_upper": (4.25, 3.95),
    "crowded": (3.95, 3.70),
    "lower_sugar": (3.75, 3.40),
    "superzoom_anomeric": (5.425, 5.385),
    "superzoom_crowded": (3.735, 3.655),
}


def corr(a: pd.Series, b: pd.Series) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[good], y[good])[0, 1])


def draw(field: str, name: str, limits: tuple[float, float], data: pd.DataFrame) -> Path:
    high, low = limits
    d = data[(data.ppm >= low) & (data.ppm <= high)].sort_values("ppm")
    ymax = float(np.nanmax(d[["experiment_norm", "spinach_norm", "gissmo_norm"]].to_numpy()))
    ymax = max(0.05, ymax * 1.12)
    fig, ax = plt.subplots(figsize=(11.5, 6.8), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(d.ppm, d.experiment_norm, color=HEDGES, lw=2.5, label="Experimental spectrum", zorder=3)
    ax.plot(d.ppm, d.spinach_norm, color=GLORY, lw=2.3, ls="--", label="Spinach simulation", zorder=4)
    ax.plot(d.ppm, d.gissmo_norm, color=OLYMPIC, lw=2.2, ls=":", label="GISSMO published simulation", zorder=2)
    ax.set_xlim(high, low)
    ax.set_ylim(-0.04 * ymax, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=14)
    ax.set_ylabel("Normalized intensity", fontsize=14)
    ax.grid(True, color=GRID, linewidth=0.9, alpha=0.85)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=12, width=1.3)
    for spine in ax.spines.values():
        spine.set_color(BLACK)
        spine.set_linewidth(1.6)

    r_se = corr(d.spinach_norm, d.experiment_norm)
    r_ge = corr(d.gissmo_norm, d.experiment_norm)
    if name == "full_sugar_region":
        title = f"Sucrose {field} MHz — full sugar region (3.00–5.80 ppm)"
        subtitle = f"Spinach vs experiment r = {r_se:.4f}   |   GISSMO vs experiment r = {r_ge:.4f}"
        ax.set_title(title + "\n" + subtitle, fontsize=17, weight="bold", pad=12)
    else:
        title = f"Sucrose {field} MHz — {name.replace('_', ' ')} ({low:.3f}–{high:.3f} ppm)"
        ax.set_title(title, fontsize=16, weight="bold", pad=12)
    leg = ax.legend(loc="upper left", fontsize=12, frameon=True, facecolor="white", edgecolor="#333333")
    leg.get_frame().set_linewidth(1.2)
    ax.text(
        0.99, 0.98,
        "Hedges: experiment\nGlory Glory dashed: Spinach\nOlympic dotted: GISSMO",
        transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
        bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.95),
    )
    fig.tight_layout()
    out = OUT / f"sucrose_{field}MHz_{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    made: list[Path] = []
    for field in FIELDS:
        path = SOURCE / f"{field}_MHz_{MODEL}_curves.csv"
        data = pd.read_csv(path)
        for name, limits in WINDOWS.items():
            made.append(draw(field, name, limits, data))
    print(f"Wrote {len(made)} Georgia-branded three-way plots to {OUT}")
    for path in made:
        print(path)


if __name__ == "__main__":
    main()
