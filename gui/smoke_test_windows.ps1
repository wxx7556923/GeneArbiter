$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest

$ProjectRoot = (Resolve-Path $PSScriptRoot).Path
$BundleDir = Join-Path $ProjectRoot "dist\GeneArbiter"
$Gui = Join-Path $BundleDir "GeneArbiter.exe"
$Worker = Join-Path $BundleDir "genearbiter-worker.exe"

if (-not (Test-Path $Gui)) {
    throw "Missing GUI executable: $Gui"
}
if (-not (Test-Path $Worker)) {
    throw "Missing worker executable: $Worker"
}

& $Worker --version | Out-Host
if ($LASTEXITCODE -ne 0) {
    throw "genearbiter-worker.exe --version failed with exit code $LASTEXITCODE"
}

$GuiProcess = Start-Process -FilePath $Gui -ArgumentList "--smoke-test" -PassThru -Wait
if ($GuiProcess.ExitCode -ne 0) {
    throw "GeneArbiter.exe --smoke-test failed with exit code $($GuiProcess.ExitCode)"
}

$SmokeOut = Join-Path ([System.IO.Path]::GetTempPath()) "GeneArbiter-build-smoke-$PID"
try {
    $Current = Join-Path $ProjectRoot "tests\fixtures\tiny_current.gff3"
    $ToolA = Join-Path $ProjectRoot "tests\fixtures\tiny_tool_a.gff3"
    $ToolB = Join-Path $ProjectRoot "tests\fixtures\tiny_tool_b.gff3"
    $Arguments = @(
        "run-files",
        "--current", $Current,
        "--candidate", "tool_a=$ToolA",
        "--candidate", "tool_b=$ToolB",
        "--out-dir", $SmokeOut,
        "--profile", "local_rule"
    )
    & $Worker @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Frozen worker tiny workflow failed with exit code $LASTEXITCODE"
    }
    foreach ($Name in @("annotation.clean.gff3", "annotation.trace.gff3", "id_mapping.tsv.gz", "run_summary.txt")) {
        $Expected = Join-Path $SmokeOut "results\$Name"
        if (-not (Test-Path $Expected)) {
            throw "Tiny workflow output is missing: $Expected"
        }
    }
} finally {
    if (Test-Path $SmokeOut) {
        Remove-Item -Recurse -Force $SmokeOut
    }
}

Write-Host "Windows GUI, worker startup, and tiny frozen workflow smoke tests passed."
