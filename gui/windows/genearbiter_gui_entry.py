"""PyInstaller entry point for the GeneArbiter desktop process."""

from genearbiter_gui.app import main


if __name__ == "__main__":
    raise SystemExit(main())
