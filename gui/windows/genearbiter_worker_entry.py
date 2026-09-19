"""PyInstaller entry point for the internal GeneArbiter console worker."""

import os

from genearbiter.cli import main


if __name__ == "__main__":
    os.environ.setdefault("GENEARBITER_IN_PROCESS_STEPS", "1")
    raise SystemExit(main())
