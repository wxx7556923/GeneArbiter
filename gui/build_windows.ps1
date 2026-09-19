param(
    [switch]$KeepBuild
)

$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = (Resolve-Path $PSScriptRoot).Path
$VenvDir = Join-Path $ProjectRoot ".venv-windows"
$VenvPython = Join-Path $VenvDir "Scripts\python.exe"
$Spec = Join-Path $ProjectRoot "windows\genearbiter_windows.spec"
$DistDir = Join-Path $ProjectRoot "dist"
$BuildDir = Join-Path $ProjectRoot "build"
$ReleaseDir = Join-Path $ProjectRoot "release"
$Archive = Join-Path $ReleaseDir "GeneArbiter-Windows-x64-0.4.0.zip"
$Checksum = "$Archive.sha256"

if (-not (Test-Path $VenvPython)) {
    if (Get-Command py -ErrorAction SilentlyContinue) {
        & py -3 -m venv $VenvDir
    } elseif (Get-Command python -ErrorAction SilentlyContinue) {
        & python -m venv $VenvDir
    } else {
        throw "Python 3.10+ was not found. Install 64-bit Python, then rerun this script."
    }
}
if (-not (Test-Path $VenvPython)) {
    throw "Virtual-environment Python was not created: $VenvPython"
}

& $VenvPython -c "import sys; print(sys.version); raise SystemExit(0 if sys.version_info >= (3, 10) else 1)"
if ($LASTEXITCODE -ne 0) {
    throw "GeneArbiter requires Python 3.10 or newer. Delete .venv-windows after installing a newer Python."
}

& $VenvPython -m pip install --upgrade pip setuptools wheel
& $VenvPython -m pip install -r (Join-Path $ProjectRoot "requirements-build.txt")
& $VenvPython -m pip install -e $ProjectRoot

Push-Location $ProjectRoot
try {
    & $VenvPython -m unittest discover -s tests -p "test_*.py" -v
    & $VenvPython -m compileall -q genearbiter genearbiter_gui
    if (-not $KeepBuild) {
        if (Test-Path $BuildDir) {
            Remove-Item -Recurse -Force $BuildDir
        }
        if (Test-Path $DistDir) {
            Remove-Item -Recurse -Force $DistDir
        }
    }
    & $VenvPython -m PyInstaller --noconfirm --clean --distpath $DistDir --workpath $BuildDir $Spec
    # Qt Virtual Keyboard is unused by this QWidget interface and has a separate GPL license.
    Remove-Item -Force -ErrorAction SilentlyContinue (Join-Path $DistDir "GeneArbiter\_internal\PySide6\Qt6VirtualKeyboard.dll")
    Remove-Item -Force -ErrorAction SilentlyContinue (Join-Path $DistDir "GeneArbiter\_internal\PySide6\plugins\platforminputcontexts\qtvirtualkeyboardplugin.dll")
    & (Join-Path $ProjectRoot "smoke_test_windows.ps1")

    $BundleDir = Join-Path $DistDir "GeneArbiter"
    Copy-Item (Join-Path $ProjectRoot "WINDOWS_USER_GUIDE.md") (Join-Path $BundleDir "README.md") -Force
    Copy-Item (Join-Path $ProjectRoot "LICENSE") (Join-Path $BundleDir "LICENSE") -Force
    Copy-Item (Join-Path $ProjectRoot "THIRD_PARTY_NOTICES.md") (Join-Path $BundleDir "THIRD_PARTY_NOTICES.md") -Force
    Copy-Item (Join-Path $ProjectRoot "LGPL-3.0.txt") (Join-Path $BundleDir "LGPL-3.0.txt") -Force
    Copy-Item (Join-Path $ProjectRoot "GPL-3.0.txt") (Join-Path $BundleDir "GPL-3.0.txt") -Force
    $PythonVersion = (& $VenvPython -c "import platform; print(platform.python_version())").Trim()
    $PackageVersions = & $VenvPython -m pip list --format=freeze
    @(
        "GeneArbiter GUI==0.4.0"
        "Python==$PythonVersion"
        $PackageVersions
    ) | Where-Object {
        $_ -notmatch "^-e " -and $_ -notmatch "@ file:"
    } | Set-Content -Encoding utf8 (Join-Path $BundleDir "build_environment.txt")
    New-Item -ItemType Directory -Force -Path $ReleaseDir | Out-Null
    if (Test-Path $Archive) {
        Remove-Item -Force $Archive
    }
    if (Test-Path $Checksum) {
        Remove-Item -Force $Checksum
    }
    Compress-Archive -Path $BundleDir -DestinationPath $Archive -CompressionLevel Optimal
    $Hash = (Get-FileHash -Algorithm SHA256 $Archive).Hash.ToLowerInvariant()
    "$Hash  $([System.IO.Path]::GetFileName($Archive))" | Set-Content -Encoding ascii $Checksum
    Write-Host ""
    Write-Host "Build complete: $Archive"
    Write-Host "Checksum: $Checksum"
} finally {
    Pop-Location
}
