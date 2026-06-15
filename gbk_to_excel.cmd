@echo off
setlocal

rem --- find a Python 3 interpreter at runtime (no hardcoded path; portable) ---
rem %~dp0 = the folder this .cmd lives in, so gbk_to_excel.py is always found beside it.
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
