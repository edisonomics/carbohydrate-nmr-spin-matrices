#!/usr/bin/env python3
"""Make four publication-quality super-zooms around the 3.55 ppm cluster."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


ROOT = Path("/Users/cece/Desktop/edison_lab/final_repo")
OUT = ROOT / "outputs/sucrose/common_reference_plots/superzoom_355"
FIELDS = (("600", "training"), ("800", "validation"),
          ("900", "training"), ("1100", "validation"))
HIGH, LOW = 3.60, 3.50
HEDGES, GLORY, OLYMPIC = "#B4BD00", "#E4002B", "#004E60"
GRID, BLACK = "#D9E2E5", "#000000"


def read_field(field: str):
    folder = ROOT / f"outputs/sucrose/{field}MHz_spinach_candidate"
    curve = next(folder.glob(f"*curves_{field}MHz.csv"))
    summary = json.loads(next(folder.glob("*summary.json")).read_text())
    offset = float(summary["ppm_offset_fitted"])
    rows = []
    with curve.open(newline="") as fh:
        for row in csv.DictReader(fh):
            p = float(row["ppm"]) - offset
            if LOW <= p <= HIGH:
                def val(key: str) -> float:
                    x = row[key].strip()
                    return float(x) if x and x.lower() != "nan" else np.nan
                rows.append((p, val("experiment"), val("candidate_spinach"), val("gissmo")))
    data = np.asarray(sorted(rows, key=lambda x: x[0]), dtype=float)
    return data, summary


def corr(a, b):
    good = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[good], b[good])[0, 1]) if good.sum() > 2 else float("nan")


def style(ax):
    ax.set_xlim(HIGH, LOW)
    ax.set_facecolor("white")
    ax.grid(True, color=GRID, linewidth=0.9, alpha=0.9)
    ax.set_axisbelow(True)
    ax.tick_params(labelsize=13, width=1.3)
    for spine in ax.spines.values():
        spine.set_color(BLACK)
        spine.set_linewidth(1.7)


def draw_one(field: str, role: str, data: np.ndarray, summary: dict) -> Path:
    ymax = max(0.08, float(np.nanmax(data[:, 1:])) * 1.14)
    fig, ax = plt.subplots(figsize=(12, 7.2), dpi=220)
    fig.patch.set_facecolor("white")
    ax.plot(data[:, 0], data[:, 1], color=HEDGES, lw=3.0, label="Experimental spectrum")
    ax.plot(data[:, 0], data[:, 2], color=GLORY, lw=2.8, ls="--", label="Candidate Spinach fit")
    ax.plot(data[:, 0], data[:, 3], color=OLYMPIC, lw=2.8, ls=":", label="GISSMO published simulation")
    style(ax)
    ax.set_ylim(-0.04 * ymax, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm), common DSS-referenced axis", fontsize=16, weight="bold")
    ax.set_ylabel("Normalized intensity", fontsize=16, weight="bold")
    r_se = corr(data[:, 2], data[:, 1])
    r_ge = corr(data[:, 3], data[:, 1])
    ax.set_title(f"Sucrose {field} MHz — 3.55 ppm super-zoom ({role})\n"
                 f"Spinach vs experiment r = {r_se:.4f}   |   GISSMO vs experiment r = {r_ge:.4f}",
                 fontsize=18, weight="bold", pad=14)
    leg = ax.legend(loc="upper left", fontsize=13, frameon=True, facecolor="white", edgecolor="#333333")
    leg.get_frame().set_linewidth(1.2)
    ax.text(0.99, 0.97, "Hedges: experiment\nGlory Glory dashed: Spinach\nOlympic dotted: GISSMO",
            transform=ax.transAxes, ha="right", va="top", fontsize=11.5,
            bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.95))
    fig.tight_layout()
    out = OUT / f"sucrose_{field}MHz_superzoom_3p55ppm.png"
    fig.savefig(out, dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    loaded = [(field, role, *read_field(field)) for field, role in FIELDS]
    paths = [draw_one(field, role, data, summary) for field, role, data, summary in loaded]

    # Also provide a compact 2x2 contact sheet for a poster overview.
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=False)
    fig.patch.set_facecolor("white")
    for ax, (field, role, data, summary) in zip(axes.flat, loaded):
        ax.plot(data[:, 0], data[:, 1], color=HEDGES, lw=2.3)
        ax.plot(data[:, 0], data[:, 2], color=GLORY, lw=2.1, ls="--")
        ax.plot(data[:, 0], data[:, 3], color=OLYMPIC, lw=2.1, ls=":")
        style(ax)
        ax.set_ylim(-0.03, max(0.08, float(np.nanmax(data[:, 1:])) * 1.12))
        ax.set_title(f"{field} MHz — {role}", fontsize=15, weight="bold")
        ax.text(0.02, 0.89,
                f"Spinach r={corr(data[:,2],data[:,1]):.4f}\nGISSMO r={corr(data[:,3],data[:,1]):.4f}",
                transform=ax.transAxes, fontsize=10.5,
                bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.93))
    for ax in axes[0, :]: ax.set_xlabel("")
    for ax in axes[:, 1]: ax.set_ylabel("")
    axes[1, 0].set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=14, weight="bold")
    axes[0, 0].set_ylabel("Normalized intensity", fontsize=14, weight="bold")
    axes[1, 0].set_ylabel("Normalized intensity", fontsize=14, weight="bold")
    handles = [
        plt.Line2D([0], [0], color=HEDGES, lw=2.5, label="Experiment"),
        plt.Line2D([0], [0], color=GLORY, lw=2.4, ls="--", label="Candidate Spinach"),
        plt.Line2D([0], [0], color=OLYMPIC, lw=2.4, ls=":", label="GISSMO"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=3,
               fontsize=12, frameon=True, facecolor="white", edgecolor="#555555")
    fig.suptitle("Sucrose 3.55 ppm super-zoom — all four fields", fontsize=21, weight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    sheet = OUT / "sucrose_all_fields_superzoom_3p55ppm.png"
    fig.savefig(sheet, dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(sheet.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {len(paths)} individual super-zooms and contact sheet to {OUT}")
    for path in paths:
        print(path)
    print(sheet)


if __name__ == "__main__":
    main()
