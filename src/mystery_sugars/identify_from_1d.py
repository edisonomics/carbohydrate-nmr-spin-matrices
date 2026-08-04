#!/usr/bin/env python3
"""Rank carbohydrate identity hypotheses from prepared multifield 1-D spectra.

This is an identity-free screening stage. It does not edit a spin matrix and
it never claims that 1-D evidence alone proves a molecular identity. The
prepared spectra are compared against reference proton-shift fingerprints
from the local BMRB/GISSMO cache and literature-informed topology metadata.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BMRB_CATALOG = Path(__file__).with_name("bmrb_candidate_catalog.json")


def _catalog_candidate(item: dict[str, Any]) -> dict[str, Any]:
    """Adapt one explore_BMRB record to the existing ranking contract."""

    centers = [float(value) for value in item.get("reference_anomeric_centers_ppm", [])]
    return {
        "id": item["candidate_id"],
        "name": item["name"],
        "class": "chebi_carbohydrate_unreviewed",
        "bmrb_entry": item.get("selected_bmrb_entry"),
        "reference_shift_file": None,
        "reference_anomeric_centers_ppm": centers,
        "expected_anomeric_clusters": len(centers) if centers else None,
        "expected_components": None,
        "forms": [],
        "topology": "unreviewed BMRB/ChEBI carbohydrate candidate",
        "notes": item.get("review_notes"),
        "review_status": item.get("review_status", "needs_review"),
        "bmrb_evidence": item.get("evidence", {}),
        "identity_warning": item.get("identity_warning", ""),
    }


def load_library(
    path: Path | None = None,
    *,
    include_bmrb_catalog: bool = False,
    bmrb_catalog_path: Path | None = None,
) -> list[dict[str, Any]]:
    path = path or Path(__file__).with_name("candidate_library.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates = list(payload.get("candidates", []))
    for candidate in candidates:
        candidate.setdefault("review_status", "approved")

    if not include_bmrb_catalog:
        return candidates
    catalog_path = bmrb_catalog_path or DEFAULT_BMRB_CATALOG
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    existing_entries = {
        candidate.get("bmrb_entry")
        for candidate in candidates
        if candidate.get("bmrb_entry")
    }
    for item in catalog.get("candidates", []):
        if item.get("selected_bmrb_entry") in existing_entries:
            continue
        adapted = _catalog_candidate(item)
        if adapted["reference_anomeric_centers_ppm"]:
            candidates.append(adapted)
    return candidates


def load_prepared(repo_root: Path, molecule: str) -> list[dict[str, Any]]:
    summary = repo_root / "outputs" / molecule / "prepared" / "preparation_summary.csv"
    if not summary.is_file():
        raise FileNotFoundError(
            f"Missing {summary}; prepare the 1-D spectra first with "
            "prepare_carbohydrate_spectra.py"
        )
    rows: list[dict[str, Any]] = []
    with summary.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            path = repo_root / "outputs" / molecule / "prepared" / row["fit_spectrum"]
            with path.open(newline="", encoding="utf-8") as spectrum:
                data = list(csv.DictReader(spectrum))
            ppm = np.array([float(item["ppm_dss"]) for item in data], dtype=float)
            value_key = "intensity_baseline_corrected"
            y = np.array([float(item[value_key]) for item in data], dtype=float)
            order = np.argsort(ppm)
            rows.append({
                "key": str(row["key"]),
                "field_mhz": float(row["field_mhz"]),
                "ppm": ppm[order],
                "intensity": y[order],
            })
    return rows


def _cluster_centers(values: Iterable[float], gap_ppm: float = 0.04) -> list[float]:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return []
    clusters: list[list[float]] = [[ordered[0]]]
    for value in ordered[1:]:
        if value - clusters[-1][-1] <= gap_ppm:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [float(np.mean(cluster)) for cluster in clusters]


def detect_anomeric_clusters(
    ppm: np.ndarray,
    intensity: np.ndarray,
    *,
    region: tuple[float, float] = (4.30, 5.80),
    water_region: tuple[float, float] = (4.65, 4.90),
    min_spacing_ppm: float = 0.008,
    cluster_gap_ppm: float = 0.04,
) -> list[float]:
    """Detect reproducible anomeric multiplet clusters in one prepared 1-D trace."""
    mask = (ppm >= region[0]) & (ppm <= region[1])
    mask &= ~((ppm >= water_region[0]) & (ppm <= water_region[1]))
    x, y = ppm[mask], np.asarray(intensity[mask], dtype=float)
    if x.size < 3 or not np.any(np.isfinite(y)):
        return []
    y = np.nan_to_num(y, nan=0.0)
    y -= np.percentile(y, 5)
    peak_height = float(np.max(y))
    if peak_height <= 0:
        return []
    threshold = max(0.015 * peak_height, float(np.percentile(y, 70)) * 0.35)
    candidates = [
        i for i in range(1, x.size - 1)
        if y[i] >= threshold and y[i] >= y[i - 1] and y[i] >= y[i + 1]
    ]
    candidates.sort(key=lambda i: float(y[i]), reverse=True)
    selected: list[int] = []
    for index in candidates:
        if all(abs(float(x[index] - x[other])) >= min_spacing_ppm for other in selected):
            selected.append(index)
    return _cluster_centers((float(x[index]) for index in selected), cluster_gap_ppm)


def _reference_anomeric_centers(repo_root: Path, candidate: dict[str, Any]) -> list[float]:
    embedded = candidate.get("reference_anomeric_centers_ppm")
    if embedded:
        return [float(value) for value in embedded]
    relative = candidate.get("reference_shift_file")
    if not relative:
        return []
    path = repo_root / str(relative)
    if not path.is_file():
        return []
    values: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("Atom_type") != "H":
                continue
            try:
                value = float(row["Val"])
            except (KeyError, TypeError, ValueError):
                continue
            if 4.30 <= value <= 5.80:
                values.append(value)
    return _cluster_centers(values, gap_ppm=0.04)


def _match_score(observed: list[float], reference: list[float], tolerance_ppm: float = 0.10) -> float | None:
    if not observed or not reference:
        return None
    distances = []
    for target in reference:
        distances.append(min(abs(target - value) for value in observed))
    scores = [math.exp(-0.5 * (distance / tolerance_ppm) ** 2) for distance in distances]
    return float(np.mean(scores))


def score_candidate(
    candidate: dict[str, Any], fields: list[dict[str, Any]], repo_root: Path
) -> dict[str, Any]:
    reference = _reference_anomeric_centers(repo_root, candidate)
    field_results: list[dict[str, Any]] = []
    scores: list[float] = []
    for field in fields:
        observed = detect_anomeric_clusters(field["ppm"], field["intensity"])
        expected_count = candidate.get("expected_anomeric_clusters")
        count_score = None
        if expected_count is not None:
            count_score = math.exp(-0.8 * abs(len(observed) - int(expected_count)))
        shift_score = _match_score(observed, reference)
        if shift_score is None and count_score is None:
            score = 0.0
        elif shift_score is None:
            score = float(count_score)
        elif count_score is None:
            score = float(shift_score)
        else:
            score = 0.70 * shift_score + 0.30 * count_score
        scores.append(score)
        field_results.append({
            "field_mhz": field["field_mhz"],
            "observed_anomeric_clusters_ppm": observed,
            "reference_anomeric_clusters_ppm": reference,
            "shift_match_score": shift_score,
            "count_score": count_score,
            "score": score,
        })
    mean_score = float(np.mean(scores)) if scores else 0.0
    reference_available = bool(reference)
    return {
        "candidate_id": candidate["id"],
        "name": candidate["name"],
        "class": candidate["class"],
        "bmrb_entry": candidate.get("bmrb_entry"),
        "forms": candidate.get("forms", []),
        "topology": candidate.get("topology"),
        "notes": candidate.get("notes"),
        "review_status": candidate.get("review_status", "approved"),
        "identity_warning": candidate.get("identity_warning", ""),
        "bmrb_evidence": candidate.get("bmrb_evidence", {}),
        "reference_available": reference_available,
        "mean_score": mean_score,
        "field_results": field_results,
    }


def rank_candidates(
    fields: list[dict[str, Any]], repo_root: Path, library: list[dict[str, Any]] | None = None
) -> list[dict[str, Any]]:
    results = [score_candidate(candidate, fields, repo_root) for candidate in (library or load_library())]
    return sorted(results, key=lambda item: item["mean_score"], reverse=True)


def build_report(repo_root: Path, molecule: str, ranked: list[dict[str, Any]]) -> dict[str, Any]:
    if not ranked:
        raise ValueError("No candidates were supplied")
    top = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    margin = float(top["mean_score"] - second["mean_score"]) if second else 1.0
    promotable = bool(
        top["reference_available"]
        and top.get("review_status", "approved") == "approved"
        and top["mean_score"] >= 0.75
        and margin >= 0.10
    )
    status = "CANDIDATE" if promotable else "REVIEW"
    return {
        "molecule": molecule,
        "stage": "identity_free_multifield_1d_screen",
        "status": status,
        "identity_claim": top["name"] if promotable else "unknown carbohydrate",
        "top_candidate": top["name"],
        "topology_hypothesis": top["topology"],
        "top_review_status": top.get("review_status", "approved"),
        "top_score": top["mean_score"],
        "runner_up": second["name"] if second else None,
        "score_margin": margin,
        "identity_confirmation_required": True,
        "confirmation_options": ["verified GISSMO/BMRB matrix", "COSY or TOCSY", "HSQC/HMBC", "matched authentic standard"],
        "ranked_candidates": ranked,
        "warning": "1-D multifield screening ranks hypotheses; it does not prove atom assignments or molecular identity.",
    }


def write_outputs(repo_root: Path, molecule: str, report: dict[str, Any]) -> tuple[Path, Path]:
    output_dir = repo_root / "outputs" / molecule / "identification"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "identification_1d_report.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    md_path = output_dir / "identification_1d_report.md"
    lines = [
        f"# 1-D mystery-sugar screen: {molecule}",
        "",
        f"Status: **{report['status']}**",
        f"Top candidate: **{report['top_candidate']}**",
        f"Top score: {report['top_score']:.3f}",
        f"Runner-up: {report['runner_up'] or 'none'}",
        f"Score margin: {report['score_margin']:.3f}",
        "",
        "This is a multifield 1-D hypothesis ranking, not an identity confirmation.",
        "",
        "| Rank | Candidate | Score | Reference shifts |",
        "|---:|---|---:|---|",
    ]
    for rank, item in enumerate(report["ranked_candidates"], 1):
        lines.append(f"| {rank} | {item['name']} | {item['mean_score']:.3f} | {'yes' if item['reference_available'] else 'no'} |")
    lines.extend(["", "Required next evidence: " + ", ".join(report["confirmation_options"]) + "."])
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--library", type=Path)
    parser.add_argument(
        "--include-bmrb-catalog",
        action="store_true",
        help="Include review-gated candidates synced from explore_BMRB",
    )
    parser.add_argument(
        "--bmrb-catalog",
        type=Path,
        default=DEFAULT_BMRB_CATALOG,
        help=f"Synced explore_BMRB catalog (default: {DEFAULT_BMRB_CATALOG})",
    )
    args = parser.parse_args()
    root = args.repo_root.resolve()
    fields = load_prepared(root, args.molecule)
    library = load_library(
        args.library,
        include_bmrb_catalog=args.include_bmrb_catalog,
        bmrb_catalog_path=args.bmrb_catalog,
    )
    ranked = rank_candidates(fields, root, library)
    report = build_report(root, args.molecule, ranked)
    json_path, md_path = write_outputs(root, args.molecule, report)
    print(f"1-D identity screen: {report['status']}")
    print(f"Top candidate: {report['top_candidate']} (score={report['top_score']:.3f})")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
