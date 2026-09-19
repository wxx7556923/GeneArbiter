from __future__ import annotations

import importlib.util
import unittest

import genearbiter
from genearbiter.cli import main


class PublicNamespaceTests(unittest.TestCase):
    def test_public_namespace_exposes_version(self) -> None:
        self.assertEqual(genearbiter.__version__, "0.4.0")

    def test_public_cli_uses_genearbiter_namespace(self) -> None:
        self.assertEqual(main.__module__, "genearbiter.cli")

    def test_legacy_package_name_is_not_shipped(self) -> None:
        legacy_name = "".join(("ai", "anno"))
        self.assertIsNone(importlib.util.find_spec(legacy_name))


if __name__ == "__main__":
    unittest.main()
