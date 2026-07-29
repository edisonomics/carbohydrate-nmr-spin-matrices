#!/usr/bin/env python3
"""Create native-resolution sucrose candidate zoom figures and montages.

The Spinach runner exports numerical curve CSVs alongside its overlay PNGs.
This utility plots those arrays directly, preserving smooth line shapes in
standalone figures.  PNG cropping is retained only as a backward-compatible
fallback for runs made before native curve export was added.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path
import zipfile

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


DEFAULT_WINDOWS = {
    "anomeric": (5.50, 5.32),
    "crowded": (3.80, 3.45),
    "fingerprint": (4.30, 3.40),
}


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _read_roles(repo: Path) -> dict[str, str]:
    path = repo / "outputs" / "sucrose" / "multifield_quality_gate.json"
    if not path.is_file():
        return {}
    data = json.loads(path.read_text())
    return {str(row["field_mhz"]): str(row["role"]) for row in data.get("fields", [])}


def _read_metrics(repo: Path) -> dict[str, tuple[float, float]]:
    path = repo / "outputs" / "sucrose" / "spinach_multifield_summary_candidate.csv"
    if not path.is_file():
        return {}
    result: dict[str, tuple[float, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            result[str(int(float(row["field_mhz"]))) ] = (
                float(row["r_spinach_vs_expt"]),
                float(row["rmse_spinach_vs_expt"]),
            )
    return result


def _source_image(repo: Path, field: str):
    path = (
        repo
        / "outputs"
        / "sucrose"
        / f"{field}MHz_spinach_candidate"
        / f"sucrose_{field}MHz_spinach_candidate_vs_expt_vs_gissmo_overlay_{field}MHz.png"
    )
    if not path.is_file():
        raise FileNotFoundError(f"Missing candidate overlay: {path}")
    image = mpimg.imread(path)
    # MATLAB export may carry transparent pixels outside the axes.  Composite
    # them onto white before placing the crop in a new figure.
    if image.ndim == 3 and image.shape[2] == 4:
        image = image.astype(float, copy=False)
        if image.max() > 1.5:  # defensive support for uint8-style readers
            image /= 255.0
        rgb = image[:, :, :3]
        alpha = np.clip(image[:, :, 3:4], 0.0, 1.0)
        image = rgb * alpha + (1.0 - alpha)
    return path, image


def _curve_source(repo: Path, field: str):
    """Locate native-resolution curves exported by the Spinach runner."""
    path = (
        repo / "outputs" / "sucrose" / f"{field}MHz_spinach_candidate"
        / f"sucrose_{field}MHz_spinach_candidate_curves_{field}MHz.csv"
    )
    if not path.is_file():
        return None
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = set(data.dtype.names or ())
    required = {"ppm", "experiment", "candidate_spinach"}
    if not required.issubset(names):
        return None
    # Older candidate runs exported the GISSMO column as all-NaN because the
    # reference curve was not passed into MATLAB.  Recover the deposited
    # GISSMO curve directly from the BMRB archive so plots remain comparable.
    if "gissmo" in names and not np.isfinite(data["gissmo"]).any():
        archive = (repo / "data" / "sucrose" / "bmrb" / "bmse000119"
                   / "bmse000119_gissmo_simulation.zip")
        member = f"bmse000119/simulation_1/B0s/sim_{field}MHz.csv"
        if archive.is_file():
            try:
                with zipfile.ZipFile(archive) as zf:
                    raw = np.genfromtxt(io.BytesIO(zf.read(member)),
                                        delimiter=",", names=True)
                xp = np.asarray(raw["ppm"], dtype=float)
                yp = np.asarray(raw["val"], dtype=float)
                order = np.argsort(xp)
                xp, yp = xp[order], yp[order]
                ppm = np.asarray(data["ppm"], dtype=float)
                gissmo = np.interp(ppm, xp, yp, left=0.0, right=0.0)
                region = (ppm >= 3.0) & (ppm <= 5.8)
                gissmo -= np.median(gissmo[region])
                anomeric = (ppm >= 5.35) & (ppm <= 5.45)
                scale = np.max(gissmo[anomeric]) if anomeric.any() else 0.0
                if np.isfinite(scale) and scale > 0:
                    data["gissmo"] = gissmo / scale
            except (KeyError, OSError, ValueError):
                pass
    return path, data


def _plot_native_curves(ax, data, window: tuple[float, float]):
    """Plot experiment, candidate, and GISSMO from one native curve table."""
    ppm = np.asarray(data["ppm"], dtype=float)
    mask = (ppm <= max(window)) & (ppm >= min(window))
    x = ppm[mask]
    if x.size == 0:
        return False

    plotted = []
    exp = np.asarray(data["experiment"], dtype=float)[mask]
    candidate = np.asarray(data["candidate_spinach"], dtype=float)[mask]
    gissmo = (np.asarray(data["gissmo"], dtype=float)[mask]
              if "gissmo" in (data.dtype.names or ()) else np.full(x.shape, np.nan))

    def add(values, color, style, label, width):
        finite = np.isfinite(values)
        if finite.any():
            line, = ax.plot(x[finite], values[finite], color=color, ls=style,
                            lw=width, label=label)
            plotted.append(line)

    add(exp, "black", "-", "Experimental spectrum", 2.0)
    add(candidate, "#d62728", "--", "Candidate Spinach fit", 2.0)
    add(gissmo, "#1769aa", ":", "Original GISSMO matrix", 2.2)
    if not plotted:
        return False

    ax.set_xlim(max(window), min(window))
    all_values = np.concatenate([
        values[np.isfinite(values)] for values in (exp, candidate, gissmo)
        if np.isfinite(values).any()
    ])
    ymin = min(-0.02, float(all_values.min()) * 1.1)
    ymax = max(0.10, float(all_values.max()) * 1.12)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=16)
    ax.set_ylabel("Normalised intensity", fontsize=16)
    ax.grid(True, color="0.88", linewidth=0.8)
    ax.tick_params(labelsize=12)
    ax.legend(loc="upper right", fontsize=12, frameon=True)
    return True


def _crop_plot(image, ppm_window: tuple[float, float]):
    """Crop a ppm window from the standard Spinach figure geometry.

    The runner uses xlim [5.8, 3.0].  The fractions below describe the
    interior axes rectangle and scale with the image dimensions, so the
    utility remains usable if export resolution changes.
    """
    height, width = image.shape[:2]
    # Interior axes rectangle measured from the standard 300-dpi export.
    x0_frac, x1_frac = 0.069, 0.985
    # Exclude the source title and the original x-axis labels.  The latter
    # become clipped when a raster crop is resized; the montage draws clean,
    # readable ppm labels below every panel instead.
    # Keep the complete white plotting background; aggressive vertical raster
    # crops can clip the source axes and create opaque black bands.  The
    # explicit legend and ppm scale below make the full-height crop readable.
    y0_frac, y1_frac = 0.180, 0.825
    x0 = width * x0_frac
    x1 = width * x1_frac
    # NMR axes run high ppm on the left to low ppm on the right.
    ppm_hi, ppm_lo = 5.8, 3.0

    def x_for_ppm(ppm: float) -> float:
        return x0 + (ppm_hi - ppm) / (ppm_hi - ppm_lo) * (x1 - x0)

    left_ppm, right_ppm = ppm_window
    left = max(0, int(round(min(x_for_ppm(left_ppm), x_for_ppm(right_ppm)) - 8)))
    right = min(width, int(round(max(x_for_ppm(left_ppm), x_for_ppm(right_ppm)) + 8)))
    top = max(0, int(round(height * (y0_frac - 0.01))))
    bottom = min(height, int(round(height * y1_frac)))
    return image[top:bottom, left:right]


def make_montage(repo: Path, name: str, window: tuple[float, float], output: Path, roles, metrics):
    fields = ["600", "800", "900", "1100"]
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), squeeze=False,
                             facecolor="white")
    fig.patch.set_facecolor("white")
    for ax, field in zip(axes.flat, fields):
        ax.set_facecolor("white")
        curves = _curve_source(repo, field)
        native = False
        if curves is not None:
            _, data = curves
            native = _plot_native_curves(ax, data, window)
        if not native:
            _, image = _source_image(repo, field)
            crop = _crop_plot(image, window)
            ax.imshow(crop, interpolation="bilinear")
            ax.set_aspect("auto")
        role = roles.get(field, "candidate")
        r, rmse = metrics.get(field, (float("nan"), float("nan")))
        ax.set_title(
            f"{field} MHz — {role} | r = {r:.4f}, RMSE = {rmse:.4f}",
            fontsize=15,
            fontweight="bold",
            pad=10,
        )
        if not native:
            ax.axis("off")
            # The source figure's legend is outside the raster crop. Repeat
            # the trace key on every fallback panel, including GISSMO.
            ax.text(
                0.02, 0.90,
                "black: experiment\nred dashed: candidate Spinach\nblue dotted: GISSMO",
                transform=ax.transAxes, fontsize=11, va="top", ha="left",
                color="black",
                bbox=dict(facecolor="white", edgecolor="0.35", alpha=0.92, pad=5),
            )
            # Draw an explicit ppm scale because the source x-axis is outside
            # the raster crop. NMR axes decrease from left to right.
            tick_ppm = np.linspace(window[0], window[1], 4)
            for i, ppm_tick in enumerate(tick_ppm):
                x_tick = i / (len(tick_ppm) - 1)
                ax.plot([x_tick, x_tick], [-0.012, 0.002], transform=ax.transAxes,
                        color="0.25", lw=1.0, clip_on=False)
                ax.text(x_tick, -0.022, f"{ppm_tick:.2f}", transform=ax.transAxes,
                        fontsize=10, ha="center", va="top", clip_on=False)
            ax.text(0.50, -0.095, r"$^1$H chemical shift (ppm)",
                    transform=ax.transAxes, fontsize=12, ha="center", va="top",
                    clip_on=False)
    fig.suptitle(
        f"Sucrose candidate Spinach fits — {name} window ({window[0]:.2f}–{window[1]:.2f} ppm)",
        fontsize=20,
        fontweight="bold",
    )
    # A figure-level legend is useful when the panels are viewed at small
    # size, while the panel legends above keep each crop self-contained.
    handles = [
        Line2D([0], [0], color="black", lw=2.2, label="Experimental spectrum"),
        Line2D([0], [0], color="red", lw=2.2, ls="--", label="Candidate Spinach fit"),
        Line2D([0], [0], color="#1769aa", lw=2.2, ls=":", label="Original GISSMO matrix"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.955),
        ncol=3,
        frameon=True,
        fontsize=13,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.90), h_pad=2.4, w_pad=1.5)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight", facecolor="white",
                edgecolor="white", transparent=False)
    plt.close(fig)
    print(f"Wrote {output}")


def make_single(repo: Path, field: str, name: str, window: tuple[float, float],
                output: Path, roles, metrics):
    """Write one field/window zoom as a standalone presentation figure."""
    role = roles.get(field, "candidate")
    r, rmse = metrics.get(field, (float("nan"), float("nan")))

    fig, ax = plt.subplots(figsize=(14, 7.5), facecolor="white")
    ax.set_facecolor("white")
    curves = _curve_source(repo, field)
    native = False
    if curves is not None:
        _, data = curves
        native = _plot_native_curves(ax, data, window)
    if native:
        native_note = "native-resolution curves"
    else:
        # Backward-compatible fallback for runs made before native curve
        # export was added.
        _, image = _source_image(repo, field)
        crop = _crop_plot(image, window)
        ax.imshow(crop, interpolation="bilinear", aspect="auto")
        ax.axis("off")
        native_note = "raster fallback"
    ax.set_title(
        f"Sucrose candidate Spinach fit — {field} MHz ({role})\n"
        f"{name} window ({window[0]:.2f}–{window[1]:.2f} ppm)  |  "
        f"r = {r:.4f}, RMSE = {rmse:.4f}",
        fontsize=20,
        fontweight="bold",
        pad=18,
    )
    if not native:
        ax.text(0.02, 0.91,
                "black: experiment\nred dashed: candidate Spinach fit\nblue dotted: GISSMO",
                transform=ax.transAxes, fontsize=14, va="top", ha="left",
                bbox=dict(facecolor="white", edgecolor="0.35", alpha=0.94, pad=7))

    if native:
        fig.text(0.99, 0.01, native_note, ha="right", va="bottom",
                 fontsize=9, color="0.35")
        fig.tight_layout(rect=(0, 0.02, 1, 0.88))
    else:
        # A clean, non-raster ppm scale below the fallback image.
        tick_ppm = np.linspace(window[0], window[1], 5)
        for i, ppm in enumerate(tick_ppm):
            x = i / (len(tick_ppm) - 1)
            ax.plot([x, x], [-0.012, 0.002], transform=ax.transAxes,
                    color="0.25", lw=1.2, clip_on=False)
            ax.text(x, -0.024, f"{ppm:.2f}", transform=ax.transAxes,
                    fontsize=13, ha="center", va="top", clip_on=False)
        ax.text(0.50, -0.105, r"$^1$H chemical shift (ppm)",
                transform=ax.transAxes, fontsize=15, ha="center", va="top",
                clip_on=False)
        handles = [
            Line2D([0], [0], color="black", lw=2.5, label="Experimental spectrum"),
            Line2D([0], [0], color="red", lw=2.5, ls="--", label="Candidate Spinach fit"),
            Line2D([0], [0], color="#1769aa", lw=2.5, ls=":", label="Original GISSMO matrix"),
        ]
        fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94),
                   ncol=3, frameon=True, fontsize=14)
        fig.tight_layout(rect=(0, 0.06, 1, 0.86))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white",
                edgecolor="white", transparent=False)
    plt.close(fig)
    print(f"Wrote {output}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=_repo_root())
    parser.add_argument(
        "--window",
        choices=sorted(DEFAULT_WINDOWS),
        action="append",
        help="Window(s) to write; defaults to all three.",
    )
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    roles = _read_roles(repo)
    metrics = _read_metrics(repo)
    selected = args.window or list(DEFAULT_WINDOWS)
    out_dir = repo / "outputs" / "sucrose" / "zoom_plots"
    for name in selected:
        lo, hi = DEFAULT_WINDOWS[name]
        output = out_dir / f"sucrose_candidate_zoom_{name}_{lo:.2f}_{hi:.2f}_ppm.png"
        make_montage(repo, name, (lo, hi), output, roles, metrics)
        for field in ("600", "800", "900", "1100"):
            single = out_dir / (
                f"sucrose_candidate_{field}MHz_zoom_{name}_"
                f"{lo:.2f}_{hi:.2f}_ppm.png"
            )
            make_single(repo, field, name, (lo, hi), single, roles, metrics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
