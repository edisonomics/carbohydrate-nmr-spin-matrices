#!/usr/bin/env python3
"""Load and validate cited, machine-readable carbohydrate-NMR knowledge."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable


DEFAULT_KNOWLEDGE_PATH = Path(__file__).with_name("carbohydrate_nmr_knowledge.json")


def load_knowledge(path: Path | None = None) -> dict[str, Any]:
    """Load the knowledge base and fail loudly on broken provenance."""

    knowledge_path = path or DEFAULT_KNOWLEDGE_PATH
    payload = json.loads(knowledge_path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("Unsupported carbohydrate-NMR knowledge schema")
    sources = payload.get("sources")
    rules = payload.get("rules")
    if not isinstance(sources, dict) or not sources:
        raise ValueError("Knowledge base has no sources")
    if not isinstance(rules, list) or not rules:
        raise ValueError("Knowledge base has no rules")

    seen: set[str] = set()
    for rule in rules:
        rule_id = str(rule.get("id", ""))
        if not rule_id or rule_id in seen:
            raise ValueError(f"Missing or duplicate literature rule id: {rule_id!r}")
        seen.add(rule_id)
        source_keys = [rule.get("source"), *rule.get("supporting_sources", [])]
        missing_sources = [key for key in source_keys if key not in sources]
        if missing_sources:
            raise ValueError(f"Rule {rule_id} has unknown sources: {missing_sources}")
        if not rule.get("locator") or not rule.get("statement"):
            raise ValueError(f"Rule {rule_id} lacks a source locator or statement")
        if not rule.get("applies_when"):
            raise ValueError(f"Rule {rule_id} lacks explicit applicability")
        if not rule.get("software_use"):
            raise ValueError(f"Rule {rule_id} lacks an allowed software use")
    return payload


def indexed_rules(knowledge: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(rule["id"]): rule for rule in knowledge["rules"]}


def source_record(knowledge: dict[str, Any], key: str) -> dict[str, Any]:
    source = dict(knowledge["sources"][key])
    source["key"] = key
    return source


def provenance_for_rule(
    knowledge: dict[str, Any], rule: dict[str, Any]
) -> dict[str, Any]:
    keys = [rule["source"], *rule.get("supporting_sources", [])]
    supporting_locators = rule.get("supporting_locators", {})
    sources = []
    for key in keys:
        source = source_record(knowledge, key)
        source["locator"] = (
            rule["locator"] if key == rule["source"] else supporting_locators.get(key)
        )
        if not source["locator"]:
            raise ValueError(f"Rule {rule['id']} lacks a locator for source {key}")
        sources.append(source)
    return {
        "rule_id": rule["id"],
        "locator": rule["locator"],
        "statement": rule["statement"],
        "sources": sources,
        "confidence": rule["confidence"],
        "software_use": rule["software_use"],
        "applies_when": rule["applies_when"],
    }


def coupling_rules_for_profile(
    knowledge: dict[str, Any], profile: str
) -> list[dict[str, Any]]:
    """Return ordered J1,2 rules for a named candidate profile."""

    rules = [
        rule
        for rule in knowledge["rules"]
        if rule.get("measurement") == "J1,2_hz"
        and profile in rule.get("profiles", [])
    ]
    order = {"alpha": 0, "beta": 1}
    return sorted(rules, key=lambda rule: order.get(str(rule.get("form")), 99))


def scoring_range(rule: dict[str, Any]) -> list[float] | None:
    """Return a quoted range or an explicitly labeled implementation range."""

    values = rule.get("expected_range") or rule.get("implementation_range")
    return [float(values[0]), float(values[1])] if values else None


def range_score(
    value: float,
    expected_range: Iterable[float],
    *,
    softness: float = 0.8,
) -> float:
    lower, upper = (float(item) for item in expected_range)
    distance = max(lower - float(value), 0.0, float(value) - upper)
    return math.exp(-0.5 * (distance / softness) ** 2)


def evidence_result(
    knowledge: dict[str, Any],
    rule: dict[str, Any],
    *,
    observed: Any,
    score: float | None,
    status: str,
    explanation: str,
) -> dict[str, Any]:
    """Build one explainable result that retains its scientific provenance."""

    if status not in {"supports", "contradicts", "not_observed", "caution"}:
        raise ValueError(f"Unsupported literature evidence status: {status}")
    return {
        **provenance_for_rule(knowledge, rule),
        "observed": observed,
        "score": score,
        "status": status,
        "explanation": explanation,
    }
