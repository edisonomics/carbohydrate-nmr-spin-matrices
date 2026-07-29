#!/usr/bin/env python3
"""Import every processed Bruker experiment in a bundle with one command."""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
from pathlib import Path


def load_importer():
    path = Path(__file__).with_name("import_bruker_dataset.py")
    spec = importlib.util.spec_from_file_location("import_bruker_dataset", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def ensure_config(repo_root: Path, molecule: str) -> None:
    config_path = repo_root / "data" / molecule / f"{molecule}_config.json"
    if config_path.is_file():
        return
    path = Path(__file__).with_name("init_carbohydrate.py")
    spec = importlib.util.spec_from_file_location("init_carbohydrate", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.create_config(repo_root, molecule)
    print(f"Created missing molecule configuration: {config_path}")


def parse_map(values: list[str], label: str) -> dict[str, str]:
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{label} must look like FIELD=VALUE, got {value!r}")
        key, item = value.split("=", 1)
        result[key.strip()] = item.strip()
    return result


def read_concentrations(root: Path) -> dict[str, str]:
    text = ""
    for name in ("README.txt", "README.md", "readme.txt", "readme.md"):
        path = root / name
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="ignore")
            break
    concentrations = {}
    headers = list(re.finditer(r"^(?P<field>\d+(?:\.\d+)?)\s*MHz\s*:", text, flags=re.I | re.M))
    for index, header in enumerate(headers):
        section_end = headers[index + 1].start() if index + 1 < len(headers) else len(text)
        section = text[header.end():section_end]
        match = re.search(r"(?P<mm>\d+(?:\.\d+)?)\s*mM", section, flags=re.I)
        if match:
            concentrations[header.group("field")] = match.group("mm")
    return concentrations


def pulse_program(experiment: Path) -> str:
    acqus = experiment / "acqus"
    if not acqus.is_file():
        return "unknown"
    text = acqus.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"##\$PULPROG=\s*<?([^>\s]+)", text)
    return match.group(1).lower() if match else "unknown"


def nucleus(experiment: Path) -> str:
    acqus = experiment / "acqus"
    if not acqus.is_file():
        return "unknown"
    text = acqus.read_text(encoding="utf-8", errors="ignore")
    match = re.search(r"##\$NUC1=\s*<?([^>\r\n]+)", text)
    return match.group(1).strip().lower() if match else "unknown"


def acqus_number(experiment: Path, key: str) -> float | None:
    acqus = experiment / "acqus"
    if not acqus.is_file():
        return None
    text = acqus.read_text(encoding="utf-8", errors="ignore")
    match = re.search(rf"##\${re.escape(key)}=\s*<?([^>\r\n]+)", text)
    if not match:
        return None
    try:
        return float(match.group(1).strip())
    except ValueError:
        return None


def field_from_metadata(experiment: Path) -> str:
    """Read the magnet field from Bruker BF1, without trusting folder names."""
    value = acqus_number(experiment, "BF1")
    if value is None:
        value = acqus_number(experiment, "SFO1")
    if value is None:
        raise ValueError(f"Could not determine field strength from {experiment / 'acqus'}")
    # Bruker records the actual resonance frequency (e.g. 499.84) in BF1;
    # project dataset keys use the nearest nominal field (500).
    return f"{round(value):g}"


def experiment_label(experiment: Path) -> str:
    """Choose a stable label from metadata or the folder only as a last resort."""
    name_file = experiment / "experiment_name.txt"
    if name_file.is_file():
        label = name_file.read_text(encoding="utf-8", errors="ignore").strip()
        if label:
            label = re.sub(r"[^A-Za-z0-9]+", "_", label).strip("_")
            return label.lower() or "experiment"
    return re.sub(r"[^A-Za-z0-9]+", "_", experiment.name).strip("_") or "experiment"


def pulse_class(experiment: Path) -> str:
    program = pulse_program(experiment)
    return "noesy" if "noesy" in program else ("1d" if program != "unknown" else "unknown")


def best_procno(experiment: Path, requested: str) -> str | None:
    """Choose a processed dataset, preferring the largest zero-filled SI.

    A Bruker experiment can contain several processed folders (for example
    pdata/1 and pdata/999). The default student workflow should choose the
    highest-resolution processed spectrum without requiring procno knowledge.
    An explicitly supplied procno remains authoritative.
    """
    pdata = experiment / "pdata"
    if requested != "1":
        candidate = pdata / requested
        return requested if (candidate / "1r").is_file() and (candidate / "procs").is_file() else None
    candidates = []
    for folder in pdata.iterdir() if pdata.is_dir() else []:
        if not folder.is_dir() or not (folder / "1r").is_file() or not (folder / "procs").is_file():
            continue
        text = (folder / "procs").read_text(encoding="utf-8", errors="ignore")
        match = re.search(r"##\$SI=\s*<?([^>\r\n]+)", text)
        try:
            si = int(float(match.group(1))) if match else 0
        except ValueError:
            si = 0
        candidates.append((si, folder.name))
    return max(candidates)[1] if candidates else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--procno", default="1")
    parser.add_argument("--sample-id", default=None)
    parser.add_argument("--tube-id", default="unassigned_tube")
    parser.add_argument("--concentration", action="append", default=[], metavar="FIELD=MM")
    parser.add_argument("--role", action="append", default=[], metavar="FIELD=ROLE")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument(
        "--nucleus", default="1H",
        help="Only import this observed nucleus (default: 1H); use all to include every nucleus",
    )
    parser.add_argument("--auto-metadata", action="store_true", help="infer concentration from README and comparable roles from concentration/sequence")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    importer = load_importer()
    ensure_config(args.repo_root, args.molecule)
    explicit_concentrations = parse_map(args.concentration, "--concentration")
    concentrations = dict(explicit_concentrations)
    roles = parse_map(args.role, "--role")
    if args.auto_metadata:
        for field, value in read_concentrations(args.source_root).items():
            concentrations.setdefault(field, value)
    sources = []
    skipped_nuclei = []
    candidates = sorted({path.parent for path in args.source_root.rglob("acqus")})
    for child in candidates:
        if not (child / "acqus").is_file():
            continue
        selected_procno = best_procno(child, args.procno)
        if selected_procno is not None:
            observed = nucleus(child)
            if args.nucleus.lower() != "all" and observed.lower() != args.nucleus.lower():
                skipped_nuclei.append((child, observed))
                continue
            sources.append((child, field_from_metadata(child), experiment_label(child), observed, selected_procno))
    if not sources:
        raise SystemExit("No Bruker experiment folders matching the requested nucleus were found")
    for child, observed in skipped_nuclei:
        print(f"Skipped {child.name}: observed nucleus {observed} (requested {args.nucleus})")

    def concentration_for(source: Path, field: str) -> str | None:
        if field in explicit_concentrations:
            return explicit_concentrations[field]
        if field in concentrations and any(char.isdigit() for char in source.name):
            # A field-only README value is safe only when the source name
            # actually identifies that field. Otherwise duplicate fields (for
            # example 2 mM and 100 mM at 500 MHz) would be mislabeled.
            if re.search(rf"(?<!\d){re.escape(field)}(?:\.\d+)?\s*mhz", source.name, flags=re.I):
                return concentrations[field]
        return None
    inferred_roles = {}
    if args.auto_metadata:
        comparable = {}
        for source, field, experiment, observed, selected_procno in sources:
            concentration = concentration_for(source, field)
            signature = (
                concentration or "unknown",
                pulse_class(source),
                observed,
            )
            comparable.setdefault(signature, []).append(field)
        for fields in comparable.values():
            if len(fields) >= 2:
                for index, field in enumerate(sorted(fields)):
                    inferred_roles[field] = "training" if index == 0 else "validation"
    for source, field, experiment, observed, selected_procno in sources:
        key = str(int(float(field))) if float(field).is_integer() else str(field)
        prior = next(
            (item for item in json.loads(
                (args.repo_root / "data" / args.molecule / f"{args.molecule}_config.json")
            .read_text(encoding="utf-8")
            ).get("datasets", []) if str(item.get("key")) == key),
            None,
        )
        # On replacement, retain curated sample metadata unless the student
        # explicitly supplied a replacement value. This makes re-importing a
        # cleaned Bruker bundle safe and reproducible.
        sample = args.sample_id or (prior or {}).get("sample_id") or f"{args.molecule}_bundle"
        concentration_text = concentration_for(source, field)
        concentration = float(concentration_text) if concentration_text is not None else None
        if concentration is None and prior:
            concentration = prior.get("concentration_mM")
        if args.auto_metadata and concentration is not None:
            sample = args.sample_id or (prior or {}).get("sample_id") or f"{args.molecule}_{concentration:g}mM"
        tube = args.tube_id
        if tube == "unassigned_tube" and prior:
            tube = prior.get("tube_id") or tube
        if field in roles:
            role = roles[field]
        elif args.replace and prior and prior.get("role"):
            role = prior["role"]
        else:
            role = inferred_roles.get(field, "unassigned")
        procno = selected_procno
        child_args = argparse.Namespace(
            molecule=args.molecule, source=source, field_mhz=float(field), experiment=experiment,
            procno=procno, sample_id=sample, tube_id=tube,
            concentration_mM=concentration, role=role, replace=args.replace, repo_root=args.repo_root,
        )
        result = importer.import_dataset(child_args)
        print(f"Imported {field} MHz from {source.name} -> {result['destination']}")
    print(f"Imported {len(sources)} dataset(s). Next: run plan_multifield_split.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
