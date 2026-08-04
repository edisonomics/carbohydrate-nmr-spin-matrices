from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "mystery_sugars" / "sync_explore_bmrb_catalog.py"


def load_module():
    spec = importlib.util.spec_from_file_location("bmrb_candidate_sync", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CandidateSyncTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.sync = load_module()

    def valid_payload(self):
        return {
            "schema_version": 1,
            "candidates": [{
                "candidate_id": "bmrb_bmse000001",
                "name": "Example",
                "review_status": "needs_review",
                "selected_bmrb_entry": "bmse000001",
                "reference_anomeric_centers_ppm": [5.1],
            }],
        }

    def test_validates_and_writes_catalog(self):
        payload = self.sync.validate_catalog(self.valid_payload())
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "catalog.json"
            self.sync.write_catalog(payload, output)
            saved = json.loads(output.read_text(encoding="utf-8"))
        self.assertEqual(saved["schema_version"], 1)

    def test_rejects_empty_catalog(self):
        with self.assertRaises(ValueError):
            self.sync.validate_catalog({"schema_version": 1, "candidates": []})


if __name__ == "__main__":
    unittest.main(verbosity=2)
