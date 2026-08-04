#!/usr/bin/env python3
"""Scientific/provenance tests for the executable literature knowledge base."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "src" / "common"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class LiteratureKnowledgeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rules = load_module("literature_rules", COMMON / "literature_rules.py")
        cls.knowledge = cls.rules.load_knowledge()
        cls.index = cls.rules.indexed_rules(cls.knowledge)

    def test_every_rule_has_resolvable_provenance(self):
        for rule in self.knowledge["rules"]:
            provenance = self.rules.provenance_for_rule(self.knowledge, rule)
            self.assertTrue(provenance["locator"])
            self.assertTrue(provenance["sources"])
            self.assertTrue(all(source.get("doi") for source in provenance["sources"]))
            self.assertTrue(all(source.get("verification") for source in provenance["sources"]))

    def test_generic_glucose_rule_is_not_applied_to_mannose(self):
        mannose = self.rules.coupling_rules_for_profile(self.knowledge, "mannose")
        ids = {rule["id"] for rule in mannose}
        self.assertEqual(ids, {
            "bubb.j12.alpha_d_mannose",
            "bubb.j12.beta_d_mannose",
        })
        self.assertNotIn("bubb.j12.beta_d_pyranose", ids)

    def test_mannose_tolerances_are_not_misrepresented_as_quoted_ranges(self):
        alpha = self.index["bubb.j12.alpha_d_mannose"]
        self.assertEqual(alpha["typical_value"], 1.6)
        self.assertNotIn("expected_range", alpha)
        self.assertIn("implementation_range", alpha)
        self.assertIn("software tolerance", alpha["implementation_note"])

    def test_duus_rules_include_applicability_and_exclusions(self):
        alpha = self.index["duus.1jch.alpha_d_pyranose"]
        self.assertEqual(alpha["applies_when"]["conformation"], "4C1")
        self.assertEqual(alpha["excludes"]["ring"], "furanose")

    def test_evidence_results_retain_source_and_observation(self):
        rule = self.index["bubb.j12.alpha_d_mannose"]
        result = self.rules.evidence_result(
            self.knowledge,
            rule,
            observed={"resolved_spacing_hz": 1.6},
            score=1.0,
            status="supports",
            explanation="Synthetic test",
        )
        self.assertEqual(result["rule_id"], rule["id"])
        self.assertEqual(result["observed"]["resolved_spacing_hz"], 1.6)
        self.assertEqual(result["sources"][0]["key"], "bubb_2003")


if __name__ == "__main__":
    unittest.main(verbosity=2)
