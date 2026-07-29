#!/usr/bin/env python3
"""Regression tests for identity-free multifield 1-D screening."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "mystery_sugars" / "identify_from_1d.py"


def load_module():
    spec = importlib.util.spec_from_file_location("mystery_sugar_identity_screen", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MysterySugarScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.screen = load_module()
        except ModuleNotFoundError as error:
            if error.name == "numpy":
                raise unittest.SkipTest("numpy is not installed in this Python environment")
            raise
        cls.np = cls.screen.np

    def synthetic_fields(self):
        ppm = self.np.linspace(3.0, 5.7, 4000)
        y = (
            self.np.exp(-((ppm - 5.182) / 0.0025) ** 2)
            + 0.9 * self.np.exp(-((ppm - 4.564) / 0.0025) ** 2)
        )
        return [
            {"field_mhz": 600.0, "ppm": ppm, "intensity": y},
            {"field_mhz": 900.0, "ppm": ppm + 0.0003, "intensity": y},
            {"field_mhz": 1100.0, "ppm": ppm - 0.0002, "intensity": y},
        ]

    def test_two_anomeric_clusters_are_detected(self):
        fields = self.synthetic_fields()
        clusters = self.screen.detect_anomeric_clusters(fields[0]["ppm"], fields[0]["intensity"])
        self.assertEqual(len(clusters), 2)
        self.assertAlmostEqual(clusters[0], 4.564, places=2)
        self.assertAlmostEqual(clusters[1], 5.182, places=2)

    def test_xylose_candidate_ranks_above_other_reference_candidates(self):
        ranked = self.screen.rank_candidates(self.synthetic_fields(), ROOT)
        self.assertEqual(ranked[0]["candidate_id"], "d_xylose")
        self.assertTrue(ranked[0]["reference_available"])
        self.assertGreater(ranked[0]["mean_score"], 0.8)

    def test_screen_never_calls_one_d_identity_confirmed(self):
        ranked = self.screen.rank_candidates(self.synthetic_fields(), ROOT)
        report = self.screen.build_report(ROOT, "mystery_sugar", ranked)
        self.assertIn(report["status"], {"CANDIDATE", "REVIEW"})
        self.assertNotEqual(report["status"], "CONFIRMED")
        self.assertTrue(report["identity_confirmation_required"])

    def test_reference_free_fructose_stays_review(self):
        fields = self.synthetic_fields()
        fructose = next(item for item in self.screen.load_library() if item["id"] == "d_fructose")
        result = self.screen.score_candidate(fructose, fields, ROOT)
        report = self.screen.build_report(ROOT, "mystery_sugar", [result])
        self.assertFalse(result["reference_available"])
        self.assertEqual(report["status"], "REVIEW")


if __name__ == "__main__":
    unittest.main(verbosity=2)
