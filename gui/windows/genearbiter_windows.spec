# -*- mode: python ; coding: utf-8 -*-
"""Build GeneArbiter.exe and its internal CLI worker in one Windows folder."""

from pathlib import Path

from PyInstaller.utils.hooks import collect_submodules

spec_dir = Path(SPECPATH).resolve()
if spec_dir.is_file():
    spec_dir = spec_dir.parent
root = spec_dir.parent
core_sources = []
for source in (root / "genearbiter").rglob("*.py"):
    destination = str(source.parent.relative_to(root))
    core_sources.append((str(source), destination))

core_hidden = collect_submodules("genearbiter")
gui_hidden = collect_submodules("genearbiter_gui")

gui = Analysis(
    [str(root / "windows" / "genearbiter_gui_entry.py")],
    pathex=[str(root)],
    binaries=[],
    datas=core_sources,
    hiddenimports=core_hidden + gui_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)
gui_pyz = PYZ(gui.pure)
gui_exe = EXE(
    gui_pyz,
    gui.scripts,
    [],
    exclude_binaries=True,
    name="GeneArbiter",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
)

worker = Analysis(
    [str(root / "windows" / "genearbiter_worker_entry.py")],
    pathex=[str(root)],
    binaries=[],
    datas=core_sources,
    hiddenimports=core_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["PySide6"],
    noarchive=False,
    optimize=0,
)
worker_pyz = PYZ(worker.pure)
worker_exe = EXE(
    worker_pyz,
    worker.scripts,
    [],
    exclude_binaries=True,
    name="genearbiter-worker",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=True,
)

bundle = COLLECT(
    gui_exe,
    worker_exe,
    gui.binaries,
    gui.datas,
    worker.binaries,
    worker.datas,
    strip=False,
    upx=False,
    name="GeneArbiter",
)
