#!/usr/bin/env python3
"""Make DSS-referenced, common-axis sucrose three-way plots.

The Spinach candidate CSVs are written on each acquisition's native model
axis.  During fitting, the experimental axis is shifted by the fitted
``ppm_offset_fitted`` so that it overlays that model axis.  Therefore the
corresponding DSS-axis coordinate is ``ppm_model - ppm_offset_fitted``.
This script applies that same transformation to all three curves without
refitting or changing the matrix.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


ROOT = Path("/Users/cece/Desktop/edison_lab/final_repo")
OUT = ROOT / "outputs/sucrose/common_reference_plots"
FIELDS = ("600", "800", "900", "1100")

# Georgia palette retained from the native Spinach overlays.
HEDGES = "#B4BD00"   # experiment
GLORY = "#E4002B"    # candidate Spinach
OLYMPIC = "#004E60"  # deposited GISSMO
BLACK = "#000000"
GRID = "#D9E2E5"

WINDOWS = {
    "full_sugar_region": (5.80, 3.00),
    "anomeric": (5.45, 5.35),
    "fructose_upper": (4.25, 3.95),
    "crowded": (3.95, 3.70),
    "lower_sugar": (3.75, 3.40),
    "superzoom_anomeric": (5.425, 5.385),
    "superzoom_crowded": (3.735, 3.655),
}


def read_curves(field: str):
    folder = ROOT / f"outputs/sucrose/{field}MHz_spinach_candidate"
    curve = next(folder.glob(f"*curves_{field}MHz.csv"))
    summary = next(folder.glob("*summary.json"))
    with summary.open() as fh:
        fit = json.load(fh)
    offset = float(fit["ppm_offset_fitted"])
    rows = []
    with curve.open(newline="") as fh:
        for row in csv.DictReader(fh):
            def value(key: str) -> float:
                raw = row[key].strip()
                return float(raw) if raw and raw.lower() != "nan" else np.nan
            rows.append((float(row["ppm"]) - offset,
                         value("experiment"),
                         value("candidate_spinach"),
                         value("gissmo")))
    rows.sort(key=lambda x: x[0])
    return np.asarray(rows, dtype=float), offset


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    good = np.isfinite(a) & np.isfinite(b)
    if good.sum() < 3:
        return float("nan")
    return float(np.corrcoef(a[good], b[good])[0, 1])


def draw(field: str, name: str, limits: tuple[float, float], rows: np.ndarray, offset: float) -> Path:
    high, low = limits
    d = rows[(rows[:, 0] >= low) & (rows[:, 0] <= high)]
    if len(d) == 0:
        raise RuntimeError(f"No points in {name} window for {field} MHz")
    ymax = float(np.nanmax(d[:, 1:]))
    ymax = max(0.05, ymax * 1.12)

    fig, ax = plt.subplots(figsize=(11.5, 6.8), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.plot(d[:, 0], d[:, 1], color=HEDGES, lw=2.6, label="Experimental spectrum", zorder=3)
    ax.plot(d[:, 0], d[:, 2], color=GLORY, lw=2.4, ls="--", label="Candidate Spinach fit", zorder=4)
    ax.plot(d[:, 0], d[:, 3], color=OLYMPIC, lw=2.3, ls=":", label="GISSMO published simulation", zorder=2)
    ax.set_xlim(high, low)
    ax.set_ylim(-0.04 * ymax, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=14)
    ax.set_ylabel("DSS-referenced normalized intensity", fontsize=14)
    ax.grid(True, color=GRID, linewidth=0.9, alpha=0.85)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=12, width=1.3)
    for spine in ax.spines.values():
        spine.set_color(BLACK)
        spine.set_linewidth(1.6)

    r_se = correlation(d[:, 2], d[:, 1])
    r_ge = correlation(d[:, 3], d[:, 1])
    if name == "full_sugar_region":
        title = f"Sucrose {field} MHz — common DSS-referenced axis\nfull sugar region (3.00–5.80 ppm)"
        subtitle = f"Spinach vs experiment r = {r_se:.4f}   |   GISSMO vs experiment r = {r_ge:.4f}"
        ax.set_title(title + "\n" + subtitle, fontsize=17, weight="bold", pad=12)
    else:
        title = f"Sucrose {field} MHz — {name.replace('_', ' ')}\ncommon DSS-referenced axis ({low:.3f}–{high:.3f} ppm)"
        ax.set_title(title, fontsize=16, weight="bold", pad=12)
    leg = ax.legend(loc="upper left", fontsize=12, frameon=True, facecolor="white", edgecolor="#333333")
    leg.get_frame().set_linewidth(1.2)
    ax.text(0.99, 0.98,
            f"Hedges: experiment\nGlory Glory dashed: Spinach\nOlympic dotted: GISSMO\naxis shift: {-offset:+.5f} ppm",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5,
            bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.95))
    fig.tight_layout()
    out = OUT / f"sucrose_{field}MHz_common_DSS_{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    made = []
    for field in FIELDS:
        rows, offset = read_curves(field)
        for name, limits in WINDOWS.items():
            made.append(draw(field, name, limits, rows, offset))
    print(f"Wrote {len(made)} common-DSS plots to {OUT}")
    for path in made:
        print(path)


if __name__ == "__main__":
    main()
