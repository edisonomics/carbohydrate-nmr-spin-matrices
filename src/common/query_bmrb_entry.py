#!/usr/bin/env python3
"""Download one BMRB NMR-STAR entry and extract proton shifts.

BMRB supplies assigned chemical shifts.  Connectivity and scalar couplings
still come from the structure/spectra and are deliberately left for the
Bubb-guided provisional-seed step.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bubb_rules import assess_config


def discover_bmrb_entry(source_root: Path) -> str:
    """Infer the most frequently referenced BMRB entry from student files."""
    counts: dict[str, int] = {}
    for path in source_root.rglob("*"):
        for entry in re.findall(r"bmse\d{6}", str(path), flags=re.I):
            normalized = entry.lower()
            counts[normalized] = counts.get(normalized, 0) + 1
    if not counts:
        raise ValueError(
            f"No BMRB ID (bmse######) was found below {source_root}; "
            "provide --entry explicitly or use a structure/GISSMO seed."
        )
    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
        choices = ", ".join(f"{entry} ({count} references)" for entry, count in ranked)
        raise ValueError(
            "More than one BMRB ID is equally represented in the source folder: "
            f"{choices}. Keep one entry's data together or provide --entry."
        )
    return ranked[0][0]


def discover_bmrb_entries_by_molecule(molecule: str) -> list[str]:
    """Search BMRB's official metabolomics molecule index by name.

    Student Bruker folders commonly contain only experiment numbers (2, 5,
    6, 7), not a BMRB accession. The standards index links molecule names to
    the corresponding bmse entries, so use it as the automatic fallback.
    """
    url = "https://bmrb.io/metabolomics/metabolomics_standards.php?dataset=metabolomics"
    request = Request(url, headers={"User-Agent": "edison-lab-carbohydrate-pipeline/1.0"})
    try:
        with urlopen(request, timeout=60) as response:
            page = response.read().decode("utf-8", errors="replace")
    except (HTTPError, URLError, TimeoutError) as error:
        raise ValueError(f"BMRB molecule search failed: {error}") from error

    target = re.sub(r"[^a-z0-9]+", "", molecule.lower())
    matches: list[tuple[int, str]] = []
    link_pattern = re.compile(
        r"<a[^>]+href=[\"'][^\"']*?id=(bmse\d{6})[^\"']*[\"'][^>]*>(.*?)</a>",
        flags=re.I | re.S,
    )
    for entry, label_html in link_pattern.findall(page):
        label = re.sub(r"<[^>]+>", " ", html.unescape(label_html))
        label = re.sub(r"\s+", " ", label).strip()
        normalized = re.sub(r"[^a-z0-9]+", "", label.lower())
        if not normalized or not target:
            continue
        if normalized == target:
            score = 0
        elif target in normalized or normalized in target:
            score = 1
        else:
            continue
        matches.append((score, entry.lower()))

    if not matches:
        raise ValueError(f"No BMRB metabolomics entry matched molecule name {molecule!r}")
    return list(dict.fromkeys(entry for _, entry in sorted(matches)))


def discover_bmrb_entry_from_manifest(repo_root: Path, molecule: str) -> str | None:
    """Reuse a previously verified accession when operating offline."""
    path = repo_root / "outputs" / molecule / "seed_selection.json"
    if not path.is_file():
        return None
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    entry = manifest.get("provenance", {}).get("bmrb_entry")
    return str(entry).lower() if entry else None


def connected_spin_blocks(matrix_path: Path) -> list[list[int]]:
    """Return one-based connected components of the nonzero J-coupling graph."""
    rows: list[list[float]] = []
    for line in matrix_path.read_text(encoding="utf-8").splitlines():
        values = line.split()
        if values:
            rows.append([float(value) for value in values])
    n = len(rows)
    if n == 0 or any(len(row) != n for row in rows):
        return []
    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if abs(rows[i][j]) > 1e-12 or abs(rows[j][i]) > 1e-12:
                neighbors[i].append(j)
                neighbors[j].append(i)
    seen: set[int] = set()
    blocks: list[list[int]] = []
    for start in range(n):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        block: list[int] = []
        while stack:
            node = stack.pop()
            block.append(node + 1)
            for neighbor in neighbors[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        blocks.append(sorted(block))
    return sorted(blocks, key=lambda block: block[0])


def find_atom_chem_shift_loop(text: str) -> tuple[list[str], list[list[str]]]:
    lines = text.splitlines()
    for start, line in enumerate(lines):
        if line.strip().lower() != "loop_":
            continue
        tags: list[str] = []
        pos = start + 1
        while pos < len(lines) and lines[pos].lstrip().startswith("_"):
            tags.append(lines[pos].strip().split()[0])
            pos += 1
        if not tags or not all(tag.lower().startswith("_atom_chem_shift.") for tag in tags):
            continue
        body: list[str] = []
        while pos < len(lines):
            stripped = lines[pos].strip().lower()
            if stripped == "stop_" or stripped.startswith("save_") or stripped.startswith("loop_"):
                break
            if stripped.startswith("#"):
                pos += 1
                continue
            body.append(lines[pos])
            pos += 1
        lexer = shlex.shlex("\n".join(body), posix=True)
        lexer.whitespace_split = True
        lexer.commenters = ""
        values = list(lexer)
        rows = [values[i : i + len(tags)] for i in range(0, len(values), len(tags))]
        rows = [row for row in rows if len(row) == len(tags)]
        return tags, rows
    raise ValueError("No _Atom_chem_shift loop found in the NMR-STAR entry")


def iter_star_loops(text: str):
    """Yield (tags, rows) for simple NMR-STAR loops used by evidence reports."""
    lines = text.splitlines()
    for start, line in enumerate(lines):
        if line.strip().lower() != "loop_":
            continue
        tags: list[str] = []
        pos = start + 1
        while pos < len(lines) and lines[pos].lstrip().startswith("_"):
            tags.append(lines[pos].strip().split()[0])
            pos += 1
        if not tags:
            continue
        body: list[str] = []
        while pos < len(lines):
            stripped = lines[pos].strip().lower()
            if stripped == "stop_" or stripped.startswith("save_") or stripped.startswith("loop_"):
                break
            if stripped.startswith("#"):
                pos += 1
                continue
            body.append(lines[pos])
            pos += 1
        lexer = shlex.shlex("\n".join(body), posix=True)
        lexer.whitespace_split = True
        lexer.commenters = ""
        values = list(lexer)
        rows = [values[i : i + len(tags)] for i in range(0, len(values), len(tags))]
        yield tags, [row for row in rows if len(row) == len(tags)]


def extract_bmrb_experiment_inventory(text: str) -> list[dict[str, str]]:
    """Extract named experiments, including COSY/TOCSY/HSQC/HMBC availability."""
    for tags, rows in iter_star_loops(text):
        if "_Experiment.Name" not in tags or "_Experiment.ID" not in tags:
            continue
        index = {tag: i for i, tag in enumerate(tags)}
        return [
            {"id": row[index["_Experiment.ID"]], "name": row[index["_Experiment.Name"]],
             "raw_data": row[index["_Experiment.Raw_data_flag"]] if "_Experiment.Raw_data_flag" in index else "?"}
            for row in rows
        ]
    return []


def extract_bmrb_peak_inventory(text: str) -> list[dict[str, object]]:
    """Summarize deposited spectral peak lists without inventing assignments."""
    return [
        {
            "spectral_peak_list_id": item["spectral_peak_list_id"],
            "experiment_name": item["experiment_name"],
            "dimensions": item["dimensions"],
            "peak_count": item["peak_count"],
            "numeric_peak_count": item["numeric_peak_count"],
        }
        for item in extract_bmrb_peak_data(text)
    ]


def _frame_scalar(frame: str, tag: str) -> str | None:
    """Read one scalar NMR-STAR tag from a saveframe."""
    match = re.search(rf"(?m)^\s*{re.escape(tag)}\s+(.+?)\s*$", frame)
    if not match:
        return None
    lexer = shlex.shlex(match.group(1), posix=True)
    lexer.whitespace_split = True
    lexer.commenters = ""
    values = list(lexer)
    return values[0] if values else None


def extract_bmrb_peak_data(text: str) -> list[dict[str, object]]:
    """Extract numeric peak coordinates and assignments from deposited peak lists.

    BMRB often lists a 2-D experiment even when it does not deposit a peak
    table.  This function therefore returns every spectral-peak saveframe,
    with an empty ``peaks`` list when no numeric ``_Peak_char`` coordinates
    are present.  Coordinates are kept in acquisition-dimension order and
    assignments are retained as labels only; no scalar couplings are inferred.
    """
    result: list[dict[str, object]] = []
    frames = re.finditer(
        r"(?ms)^save_spectral_peak_[^\n]*\n(.*?)(?=^save_|\Z)", text
    )
    for frame_match in frames:
        frame = frame_match.group(0)
        peak_list_id = _frame_scalar(frame, "_Spectral_peak_list.ID")
        experiment_name = _frame_scalar(frame, "_Spectral_peak_list.Experiment_name")
        dimensions = _frame_scalar(frame, "_Spectral_peak_list.Number_of_spectral_dimensions")
        if not peak_list_id or not experiment_name:
            continue

        dimension_metadata: dict[str, dict[str, object]] = {}
        coordinates: dict[str, dict[str, dict[str, object]]] = {}
        coupling_patterns: dict[str, dict[str, str]] = {}
        assignments: dict[str, dict[str, list[str]]] = {}
        for tags, rows in iter_star_loops(frame):
            index = {tag: i for i, tag in enumerate(tags)}
            if "_Spectral_dim.ID" in index:
                for row in rows:
                    dim_id = row[index["_Spectral_dim.ID"]]
                    metadata: dict[str, object] = {"id": dim_id}
                    for tag, key in (
                        ("_Spectral_dim.Atom_type", "atom_type"),
                        ("_Spectral_dim.Atom_isotope_number", "isotope"),
                        ("_Spectral_dim.Spectral_region", "region"),
                    ):
                        if tag in index:
                            metadata[key] = row[index[tag]]
                    dimension_metadata[dim_id] = metadata
            if "_Peak_char.Peak_ID" in index and "_Peak_char.Spectral_dim_ID" in index:
                value_tag = "_Peak_char.Chem_shift_val"
                pattern_tag = "_Peak_char.Coupling_pattern"
                if value_tag not in index:
                    continue
                for row in rows:
                    value_text = row[index[value_tag]]
                    if value_text in {"", ".", "?"}:
                        continue
                    try:
                        value = float(value_text)
                    except ValueError:
                        continue
                    peak_id = row[index["_Peak_char.Peak_ID"]]
                    dim_id = row[index["_Peak_char.Spectral_dim_ID"]]
                    coordinates.setdefault(peak_id, {})[dim_id] = {
                        "dimension_id": dim_id,
                        "value": value,
                    }
                    if pattern_tag in index and row[index[pattern_tag]] not in {"", ".", "?"}:
                        coupling_patterns.setdefault(peak_id, {})[dim_id] = row[index[pattern_tag]]
            if "_Assigned_peak_chem_shift.Peak_ID" in index and "_Assigned_peak_chem_shift.Spectral_dim_ID" in index:
                atom_tag = "_Assigned_peak_chem_shift.Atom_ID"
                if atom_tag not in index:
                    continue
                for row in rows:
                    atom_id = row[index[atom_tag]]
                    if atom_id in {"", ".", "?"}:
                        continue
                    peak_id = row[index["_Assigned_peak_chem_shift.Peak_ID"]]
                    dim_id = row[index["_Assigned_peak_chem_shift.Spectral_dim_ID"]]
                    assignments.setdefault(peak_id, {}).setdefault(dim_id, []).append(atom_id)

        peaks: list[dict[str, object]] = []
        def peak_sort_key(value: str) -> tuple[int, str]:
            try:
                return (0, f"{int(value):012d}")
            except ValueError:
                return (1, value)

        for peak_id in sorted(coordinates, key=peak_sort_key):
            dims = sorted(
                coordinates[peak_id],
                key=lambda value: (0, int(value)) if value.isdigit() else (1, value),
            )
            peak_coordinates: list[dict[str, object]] = []
            for dim_id in dims:
                coordinate = dict(coordinates[peak_id][dim_id])
                coordinate.update(dimension_metadata.get(dim_id, {}))
                if peak_id in coupling_patterns and dim_id in coupling_patterns[peak_id]:
                    coordinate["coupling_pattern"] = coupling_patterns[peak_id][dim_id]
                peak_coordinates.append(coordinate)
            peak: dict[str, object] = {
                "id": peak_id,
                "coordinates": peak_coordinates,
                "coordinate_values": [item["value"] for item in peak_coordinates],
            }
            if peak_id in assignments:
                peak["assigned_atoms"] = [
                    {"dimension_id": dim_id, "atom_ids": atom_ids}
                    for dim_id, atom_ids in sorted(assignments[peak_id].items())
                ]
            peaks.append(peak)

        result.append({
            "spectral_peak_list_id": peak_list_id,
            "experiment_name": experiment_name,
            "dimensions": dimensions or "?",
            "peak_count": len(peaks),
            "numeric_peak_count": sum(bool(item.get("coordinates")) for item in peaks),
            "dimension_metadata": list(dimension_metadata.values()),
            "peaks": peaks,
        })
    return result


def clean_tag(tag: str) -> str:
    return tag.rsplit(".", 1)[-1]


def _cells(row_html: str) -> list[str]:
    values = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row_html, flags=re.I | re.S)
    return [re.sub(r"<[^>]+>", "", html.unescape(value)).strip() for value in values]


def extract_gissmo_matrix(page_html: str) -> tuple[list[str], list[list[float]]] | None:
    """Extract the displayed upper-triangular GISSMO matrix from the entry page."""
    marker = re.search(r"Spin System Matrix", page_html, flags=re.I)
    if not marker:
        return None
    table_match = re.search(r"<table[^>]*>(.*?)</table>", page_html[marker.end() :], flags=re.I | re.S)
    if not table_match:
        return None
    rows = [_cells(row) for row in re.findall(r"<tr[^>]*>(.*?)</tr>", table_match.group(1), flags=re.I | re.S)]
    rows = [row for row in rows if row]
    if len(rows) < 2:
        return None
    atom_ids = rows[0][1:]
    matrix: list[list[float]] = []
    for row in rows[1:]:
        if len(row) != len(atom_ids) + 1:
            continue
        try:
            matrix.append([float(value) for value in row[1:]])
        except ValueError:
            continue
    if len(matrix) != len(atom_ids) or not atom_ids:
        return None
    return atom_ids, matrix


def download_gissmo(entry: str | list[str], root: Path) -> dict:
    """Download GISSMO project files and the displayed numeric matrix if present."""
    entries = [entry] if isinstance(entry, str) else entry
    entries = list(dict.fromkeys(entries))
    result: dict = {"available": False, "checked_entries": entries}
    for candidate in entries:
        base = f"https://gissmo.bmrb.io/entry/{candidate}/simulation_1"
        for suffix, filename in (("/spin_simulation.xml", f"{candidate}_spin_simulation.xml"),
                                 (f"/{candidate}_simulation_1_nmredata.zip", f"{candidate}_simulation_1_nmredata.zip"),
                                 ("/zip", f"{candidate}_gissmo_simulation.zip")):
            url = base + suffix
            try:
                request = Request(url, headers={"User-Agent": "edison-lab-carbohydrate-pipeline/1.0"})
                with urlopen(request, timeout=60) as response:
                    payload = response.read()
                # GISSMO sometimes returns HTTP 200 with a short text error
                # such as "No such entry exists." rather than a real artifact.
                # Never record that as an available matrix/simulation file.
                stripped = payload.lstrip()
                if stripped.lower().startswith(b"no such entry"):
                    continue
                if filename.endswith(".xml") and not stripped.startswith(b"<"):
                    continue
                if filename.endswith(".zip") and not payload.startswith(b"PK"):
                    continue
                (root / filename).write_bytes(payload)
                result[filename] = str(root / filename)
                result["available"] = True
            except (HTTPError, URLError, TimeoutError):
                continue
        try:
            request = Request(base, headers={"User-Agent": "edison-lab-carbohydrate-pipeline/1.0"})
            with urlopen(request, timeout=60) as response:
                page = response.read().decode("utf-8", errors="replace")
            parsed = extract_gissmo_matrix(page)
            if parsed is not None:
                atom_ids, matrix = parsed
                matrix_path = root / "gissmo_spin_matrix.txt"
                matrix_path.write_text("\n".join(" ".join(f"{value:.9f}" for value in row) for row in matrix) + "\n", encoding="utf-8")
                (root / "gissmo_atom_ids.json").write_text(json.dumps(atom_ids, indent=2) + "\n", encoding="utf-8")
                result["matrix_file"] = str(matrix_path)
                result["atom_ids"] = atom_ids
                result["entry_url"] = base
                result["gissmo_entry"] = candidate
                result["available"] = True
                return result
        except (HTTPError, URLError, TimeoutError):
            continue
    return result


def fetch_bmrb_entry(entry: str) -> tuple[bytes, str]:
    """Fetch a metabolite entry through REST, then the official FTP archive."""
    urls = [
        f"https://bmrb.io/rest/bmrb/{entry}/nmr-star3",
        f"https://bmrb.io/ftp/pub/bmrb/metabolomics/entry_directories/{entry}/{entry}.str",
    ]
    last_error: Exception | None = None
    for url in urls:
        try:
            request = Request(url, headers={"User-Agent": "edison-lab-carbohydrate-pipeline/1.0"})
            with urlopen(request, timeout=60) as response:
                return response.read(), url
        except (HTTPError, URLError, TimeoutError) as error:
            last_error = error
    raise RuntimeError(f"BMRB entry {entry} could not be downloaded from REST or metabolomics archive: {last_error}")


def related_bmrb_entries(entry: str) -> list[str]:
    """Find related metabolomics IDs listed on the BMRB compound page."""
    url = f"https://bmrb.io/metabolomics/mol_summary/show_data.php?id={entry}&whichTab=1"
    try:
        request = Request(url, headers={"User-Agent": "edison-lab-carbohydrate-pipeline/1.0"})
        with urlopen(request, timeout=60) as response:
            page = response.read().decode("utf-8", errors="replace")
        ids = re.findall(r"bmse\d{6}", page, flags=re.I)
        return list(dict.fromkeys([entry] + ids))
    except (HTTPError, URLError, TimeoutError):
        return [entry]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--entry", default="auto",
        help="BMRB ID, e.g. bmse000119; use auto to infer it from files or molecule name",
    )
    parser.add_argument("--molecule", required=True)
    parser.add_argument(
        "--source-root", type=Path,
        help="Optional student data folder used to detect a BMRB ID",
    )
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()

    entry = args.entry.strip()
    if entry.lower() == "auto":
        source_error: ValueError | None = None
        if args.source_root is not None:
            try:
                entry = discover_bmrb_entry(args.source_root.expanduser().resolve())
                print(f"Auto-selected BMRB entry {entry} from {args.source_root}")
            except ValueError as error:
                source_error = error
        if not entry or entry.lower() == "auto":
            manifest_entry = discover_bmrb_entry_from_manifest(args.repo_root, args.molecule)
            if manifest_entry:
                entry = manifest_entry
                print(f"Auto-selected previously verified BMRB entry {entry} from the seed manifest")
            else:
                try:
                    candidates = discover_bmrb_entries_by_molecule(args.molecule)
                except ValueError as error:
                    detail = f"; file search also failed: {source_error}" if source_error else ""
                    raise SystemExit(str(error) + detail) from error
                entry = candidates[0]
                if len(candidates) > 1:
                    print(f"BMRB molecule search candidates: {', '.join(candidates)}")
                print(f"Auto-selected BMRB entry {entry} by molecule name: {args.molecule}")
    try:
        raw, bmrb_url = fetch_bmrb_entry(entry)
    except RuntimeError as error:
        cached = args.repo_root / "data" / args.molecule / "bmrb" / entry / f"{entry}.str"
        if cached.is_file():
            raw = cached.read_bytes()
            bmrb_url = f"cached:{cached}"
            print(f"BMRB network unavailable; using cached entry {cached}")
        else:
            raise SystemExit(str(error)) from error
    if not raw.strip():
        raise RuntimeError(f"BMRB returned an empty entry: {bmrb_url}")

    root = args.output_dir or args.repo_root / "data" / args.molecule / "bmrb" / entry
    root.mkdir(parents=True, exist_ok=True)
    star_path = root / f"{entry}.str"
    star_path.write_bytes(raw)
    tags, rows = find_atom_chem_shift_loop(raw.decode("utf-8", errors="replace"))
    columns = [clean_tag(tag) for tag in tags]
    records = [dict(zip(columns, row)) for row in rows]
    proton_rows = [row for row in records if row.get("Atom_type", "").upper() == "H"]
    csv_path = root / "chemical_shifts.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(proton_rows)

    atoms = []
    for row in proton_rows:
        value = row.get("Val", "")
        atom_id = row.get("Atom_ID", "")
        if value in {"", ".", "?"} or not atom_id:
            continue
        try:
            shift = float(value)
        except ValueError:
            continue
        atoms.append({
            "id": atom_id,
            "shift_ppm": shift,
            "comp_id": row.get("Comp_ID"),
            "comp_index_id": row.get("Comp_index_ID"),
            "shift_sigma_ppm": None,
            "source": "BMRB",
        })
    star_text = raw.decode("utf-8", errors="replace")
    experiments = extract_bmrb_experiment_inventory(star_text)
    peak_data = extract_bmrb_peak_data(star_text)
    peak_lists = [
        {
            "spectral_peak_list_id": item["spectral_peak_list_id"],
            "experiment_name": item["experiment_name"],
            "dimensions": item["dimensions"],
            "peak_count": item["peak_count"],
            "numeric_peak_count": item["numeric_peak_count"],
        }
        for item in peak_data
    ]
    observations = {
        "molecule": args.molecule,
        "atoms": atoms,
        "couplings": [],
        "experiments": experiments,
        "peak_lists": peak_lists,
        "peak_data": peak_data,
        "evidence": [{
            "type": "BMRB_experiment_inventory",
            "experiment_count": len(experiments),
            "peak_list_count": len(peak_lists),
            "numeric_peak_list_count": sum(bool(item["numeric_peak_count"]) for item in peak_data),
            "numeric_2d_peak_list_count": sum(
                str(item.get("dimensions")) == "2" and bool(item["numeric_peak_count"])
                for item in peak_data
            ),
            "note": "Numeric peak coordinates and deposited assignments are preserved; scalar couplings are never inferred from peak positions.",
        }],
        "notes": "BMRB shifts and deposited peak coordinates imported; add Bubb-guided topology and measured COSY/TOCSY/HSQC couplings before matrix refinement.",
    }
    observations_path = root / "spectral_observations.json"
    observations_path.write_text(json.dumps(observations, indent=2) + "\n", encoding="utf-8")
    gissmo = download_gissmo(related_bmrb_entries(entry), root)
    cached_matrix = root / "gissmo_spin_matrix.txt"
    if not gissmo.get("matrix_file") and cached_matrix.is_file():
        cached_ids = root / "gissmo_atom_ids.json"
        gissmo = {
            "available": True,
            "matrix_file": str(cached_matrix),
            "atom_ids": json.loads(cached_ids.read_text(encoding="utf-8")) if cached_ids.is_file() else [],
            "entry_url": f"https://gissmo.bmrb.io/entry/{entry}/simulation_1",
            "gissmo_entry": entry,
            "cached": True,
        }
        print(f"GISSMO network unavailable; using cached matrix {cached_matrix}")
    next_step = (
        "Run multifield transfer and identifiability validation on the downloaded GISSMO matrix."
        if gissmo.get("matrix_file")
        else "Add Bubb-guided connectivity and measured COSY/TOCSY/HSQC couplings, then run build_provisional_seed.py."
    )
    provenance = {
        "source": "BMRB",
        "entry": entry,
        "url": bmrb_url,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "raw_nmr_star": str(star_path),
        "chemical_shifts_csv": str(csv_path),
        "proton_shift_count": len(atoms),
        "experiment_count": len(experiments),
        "peak_list_count": len(peak_lists),
        "numeric_peak_list_count": sum(bool(item["numeric_peak_count"]) for item in peak_data),
        "numeric_2d_peak_list_count": sum(
            str(item.get("dimensions")) == "2" and bool(item["numeric_peak_count"])
            for item in peak_data
        ),
        "two_d_peak_count": sum(
            int(item["peak_count"])
            for item in peak_data
            if str(item.get("dimensions")) == "2"
        ),
        "two_d_experiments": [
            item["name"] for item in experiments
            if any(token in item["name"].upper() for token in ("COSY", "TOCSY", "HSQC", "HMBC"))
        ],
        "gissmo": gissmo,
        "next_step": next_step,
    }
    (root / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    if gissmo.get("matrix_file"):
        matrix_path = Path(gissmo["matrix_file"])
        try:
            relative_matrix = str(matrix_path.resolve().relative_to(args.repo_root.resolve()))
        except ValueError:
            relative_matrix = str(matrix_path.resolve())
        seed_manifest = {
            "molecule": args.molecule,
            "requested_source": "bmrb",
            "selected_source": "gissmo_bmrb_endpoint",
            "status": "READY",
            "confidence": "verified",
            "matrix_file": relative_matrix,
            "atom_ids": gissmo.get("atom_ids", []),
            "blocks": connected_spin_blocks(Path(gissmo["matrix_file"])),
            "provenance_required": True,
            "provenance": {"bmrb_entry": entry, "gissmo_entry_url": gissmo["entry_url"]},
            "next_step": "Run multifield transfer and identifiability validation before publication.",
            "created_utc": datetime.now(timezone.utc).isoformat(),
        }
        seed_path = args.repo_root / "outputs" / args.molecule / "seed_selection.json"
        seed_path.parent.mkdir(parents=True, exist_ok=True)
        seed_path.write_text(json.dumps(seed_manifest, indent=2) + "\n", encoding="utf-8")
        config_path = args.repo_root / "data" / args.molecule / f"{args.molecule}_config.json"
        if config_path.is_file():
            config = json.loads(config_path.read_text(encoding="utf-8"))
            config["matrix_file"] = relative_matrix
            config["atom_ids"] = gissmo.get("atom_ids", config.get("atom_ids", []))
            derived_blocks = connected_spin_blocks(matrix_path)
            if derived_blocks:
                config["blocks"] = derived_blocks
            config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            print(f"Updated molecule configuration: {config_path}")
        print(f"GISSMO matrix found and selected: {matrix_path}")
        print(f"Wrote seed manifest: {seed_path}")
    else:
        print("No GISSMO matrix was found for this BMRB entry; using the shift observations as the provisional-seed input.")
        print("Next: review COSY/TOCSY/HSQC assignments and run the independent J-measurement stage before matrix refinement.")

    # Make the chemistry interpretation explicit for students and reviewers.
    # Bubb supplies the guardrails; BMRB/GISSMO supplies the molecule-specific
    # shifts and matrix.  This report is advisory and never invents couplings.
    config_path = args.repo_root / "data" / args.molecule / f"{args.molecule}_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else {"name": args.molecule}
    chemistry_report = assess_config(
        args.molecule,
        config,
        bmrb={
            "proton_shift_count": len(atoms),
            "gissmo_matrix_file": gissmo.get("matrix_file"),
            "two_d_experiments": provenance.get("two_d_experiments", []),
            "two_d_peak_count": provenance.get("two_d_peak_count", 0),
            "entry": entry,
            "source": "BMRB",
        },
    )
    chemistry_report["bmrb_entry"] = entry
    chemistry_report["bubb_reference"] = chemistry_report["profile"]["reference"]
    chemistry_path = args.repo_root / "outputs" / args.molecule / "chemistry_evidence.json"
    chemistry_path.parent.mkdir(parents=True, exist_ok=True)
    chemistry_path.write_text(json.dumps(chemistry_report, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote Bubb/BMRB chemistry report: {chemistry_path}")
    print(f"Downloaded {entry} from {bmrb_url}")
    print(f"Extracted {len(atoms)} proton shifts")
    print(f"Wrote {star_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {observations_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
