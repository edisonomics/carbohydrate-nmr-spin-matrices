#!/usr/bin/env python3
"""Sync the review-gated candidate catalog exported by explore_BMRB."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Any
from urllib.request import Request, urlopen


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOCAL_SOURCE = REPO_ROOT.parent / "explore_BMRB" / "data" / "bmrb_mystery_sugar_candidates.json"
DEFAULT_REMOTE_SOURCE = (
    "https://piqueen314.github.io/explore_BMRB/"
    "data/downloads/bmrb_mystery_sugar_candidates.json"
)
DEFAULT_OUTPUT = Path(__file__).with_name("bmrb_candidate_catalog.json")


def validate_catalog(payload: Any) -> dict[str, Any]:
    """Validate the small contract consumed by the mystery-sugar ranker."""

    if not isinstance(payload, dict) or payload.get("schema_version") != 1:
        raise ValueError("Expected explore_BMRB candidate catalog schema version 1")
    candidates = payload.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("The explore_BMRB candidate catalog contains no candidates")
    required = {
        "candidate_id",
        "name",
        "review_status",
        "selected_bmrb_entry",
        "reference_anomeric_centers_ppm",
    }
    for index, candidate in enumerate(candidates):
        if not isinstance(candidate, dict):
            raise ValueError(f"Candidate {index} is not an object")
        missing = required - candidate.keys()
        if missing:
            raise ValueError(
                f"Candidate {index} is missing fields: {', '.join(sorted(missing))}"
            )
    return payload


def load_catalog(source: Path | None = None, url: str = DEFAULT_REMOTE_SOURCE) -> tuple[dict[str, Any], str]:
    """Load a sibling-repository export when present, otherwise the live copy."""

    local_source = source or DEFAULT_LOCAL_SOURCE
    if local_source.is_file():
        return validate_catalog(json.loads(local_source.read_text(encoding="utf-8"))), str(local_source)

    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "carbohydrate-nmr-spin-matrices/1.0",
        },
    )
    with urlopen(request, timeout=120.0) as response:
        return validate_catalog(json.load(response)), url


def write_catalog(payload: dict[str, Any], output: Path = DEFAULT_OUTPUT) -> None:
    """Write the synced catalog atomically."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary: Path | None = None
    try:
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=output.parent,
            prefix=f".{output.name}.",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        help="Local explore_BMRB JSON export (defaults to the sibling repository)",
    )
    parser.add_argument("--url", default=DEFAULT_REMOTE_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    payload, source = load_catalog(args.source, args.url)
    write_catalog(payload, args.output)
    summary = payload.get("summary", {})
    print(
        f"Synced {len(payload['candidates'])} structures from {source} to "
        f"{args.output}; {summary.get('candidates_with_anomeric_references', 0)} "
        "have anomeric-region references"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
