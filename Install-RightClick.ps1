<#
    Install-RightClick.ps1

    Adds a right-click "Extract annotations to Excel" option for GenBank files
    (.gbk / .gbff / .gb / .genbank) and for folders containing them.

    - Per-user only (writes to HKCU and your SendTo folder). NO admin rights needed.
    - Fully reversible with Uninstall-RightClick.ps1.

    How to run:
        Right-click this file  ->  "Run with PowerShell"
        ...or in a terminal:   powershell -ExecutionPolicy Bypass -File .\Install-RightClick.ps1

    Where the option appears:
      * Files/folders: right-click -> (on Windows 11) "Show more options" ->
                       "Extract annotations to Excel".
      * Also added under right-click -> "Send to" -> "Extract GBK to Excel"
        (this one shows in the normal Windows 11 menu too).
#>

$ErrorActionPreference = "Stop"
$scriptDir = $PSScriptRoot
$pyScript  = Join-Path $scriptDir "gbk_to_excel.py"

if (-not (Test-Path $pyScript)) {
    throw "gbk_to_excel.py was not found next to this installer ($scriptDir)."
}

# --- locate Python -------------------------------------------------------
$python = (Get-Command python -ErrorAction SilentlyContinue).Source
if (-not $python) { $python = (Get-Command python3 -ErrorAction SilentlyContinue).Source }
if (-not $python) {
    $pyLauncher = (Get-Command py -ErrorAction SilentlyContinue).Source
    if ($pyLauncher) { $python = & $pyLauncher -3 -c "import sys; print(sys.executable)" }
}
if (-not $python) {
    throw "Python was not found on PATH. Install Python 3, or edit this script with its path."
}
Write-Host "Using Python: $python"

# --- ensure openpyxl is available ---------------------------------------
& $python -c "import openpyxl" 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "Installing the 'openpyxl' package (needed to write Excel files)..."
    & $python -m pip install --quiet openpyxl
}

# --- write the wrapper .cmd ---------------------------------------------
# The wrapper finds Python itself at runtime (py launcher / PATH) rather than
# baking in this machine's interpreter path, so the file stays portable: it can
# be copied to any machine and works without re-running this installer.
# %~dp0 = folder this .cmd lives in, so it always finds gbk_to_excel.py beside it.
$cmdPath = Join-Path $scriptDir "gbk_to_excel.cmd"
$cmdBody = @'
@echo off
setlocal

rem --- find a Python 3 interpreter at runtime (no hardcoded path; portable) ---
set "PY="
where py >nul 2>nul && set "PY=py -3"
if not defined PY (
  where python >nul 2>nul && set "PY=python"
)
if not defined PY (
  echo.
  echo *** Python 3 was not found on PATH. ***
  echo     Install it from https://www.python.org/downloads/ ^(tick "Add python.exe to PATH"^),
  echo     then run Install-RightClick.ps1 again.
  pause
  exit /b 9009
)

%PY% "%~dp0gbk_to_excel.py" --open %*
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" (
  echo.
  echo *** gbk_to_excel failed ^(exit code %RC%^). Check the messages above. ***
  pause
)
endlocal
'@
Set-Content -Path $cmdPath -Value $cmdBody -Encoding ascii
Write-Host "Wrote wrapper: $cmdPath"

$command = "`"$cmdPath`" `"%1`""

# --- register context-menu verbs (HKCU) --------------------------------
function Add-Verb($keyPath, $label) {
    New-Item -Path "$keyPath\command" -Force | Out-Null
    Set-Item  -Path $keyPath -Value $label
    New-ItemProperty -Path $keyPath -Name "Icon" -Value $python -PropertyType String -Force | Out-Null
    Set-Item  -Path "$keyPath\command" -Value $command
}

foreach ($ext in @(".gbk", ".gbff", ".gb", ".genbank")) {
    Add-Verb "HKCU:\Software\Classes\SystemFileAssociations\$ext\shell\GbkToExcel" "Extract annotations to Excel"
}
Add-Verb "HKCU:\Software\Classes\Directory\shell\GbkToExcel" "Extract GBK annotations to Excel"
Write-Host "Registered right-click menu for .gbk/.gbff/.gb/.genbank files and folders."

# --- Send To shortcut (visible in the normal Windows 11 menu) -----------
$sendTo = Join-Path $env:APPDATA "Microsoft\Windows\SendTo"
$ws  = New-Object -ComObject WScript.Shell
$lnk = $ws.CreateShortcut((Join-Path $sendTo "Extract GBK to Excel.lnk"))
$lnk.TargetPath       = $cmdPath
$lnk.WorkingDirectory = $scriptDir
$lnk.IconLocation     = "$python,0"
$lnk.Description       = "Extract GenBank annotations into an Excel workbook"
$lnk.Save()
Write-Host "Added 'Send to' -> 'Extract GBK to Excel'."

Write-Host ""
Write-Host "DONE. Try it now:" -ForegroundColor Green
Write-Host "  Right-click any .gbk file  ->  Send to  ->  Extract GBK to Excel"
Write-Host "  (or: right-click -> Show more options -> Extract annotations to Excel)"
Write-Host ""
Write-Host "To remove it later, run Uninstall-RightClick.ps1."
