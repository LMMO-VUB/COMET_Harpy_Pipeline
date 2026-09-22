@echo off
setlocal enabledelayedexpansion
title COMET Pipeline - New Run Setup

echo ============================================================
echo   COMET Spatial Pipeline - Set Up a New Run Folder
echo ============================================================
echo.
echo This builds a run folder with config.yaml already filled in,
echo so you don't have to write one by hand.
echo.
echo Drag your HALO or Horizon export CSV onto this window and press
echo Enter, or type/paste its full path below.
echo.

set "SRC_CSV=%~1"
if "%SRC_CSV%"=="" (
    set /p SRC_CSV=Path to your HALO/Horizon export CSV:
)
set SRC_CSV=%SRC_CSV:"=%

if "%SRC_CSV%"=="" (
    echo.
    echo ERROR: No file given.
    echo.
    pause
    exit /b 1
)

if not exist "%SRC_CSV%" (
    echo.
    echo ERROR: Could not find the file:
    echo   "%SRC_CSV%"
    echo.
    pause
    exit /b 1
)

set "CSV_NAME="
set "DEFAULT_PROJECT="
set "SRC_DIR="
for %%F in ("%SRC_CSV%") do (
    set "CSV_NAME=%%~nxF"
    set "DEFAULT_PROJECT=%%~nF"
    set "SRC_DIR=%%~dpF"
)

echo.
set /p PROJECT=Project name [%DEFAULT_PROJECT%]:
if "%PROJECT%"=="" set "PROJECT=%DEFAULT_PROJECT%"
rem Strip apostrophes - they would break the generated PowerShell script below.
set "PROJECT=%PROJECT:'=%"

set "RUNFOLDER=%SRC_DIR%%PROJECT%_run"

if exist "%RUNFOLDER%\" (
    echo.
    echo ERROR: "%RUNFOLDER%" already exists.
    echo Pick a different project name, or remove/rename that folder
    echo and run this again.
    echo.
    pause
    exit /b 1
)

set "TEMPLATE=%~dp0config.yaml.template"
if not exist "%TEMPLATE%" (
    echo.
    echo ERROR: config.yaml.template not found next to this script in
    echo   "%~dp0"
    echo Ask whoever maintains the pipeline to add it to the
    echo pipeline_code folder.
    echo.
    pause
    exit /b 1
)

echo.
echo Creating run folder ...
mkdir "%RUNFOLDER%"

echo Copying your CSV into the run folder ^(this can take a moment
echo for large exports^) ...
copy "%SRC_CSV%" "%RUNFOLDER%\%CSV_NAME%" >nul
if errorlevel 1 (
    echo.
    echo ERROR: Failed to copy the CSV into "%RUNFOLDER%".
    echo.
    pause
    exit /b 1
)

echo Writing config.yaml ...
set "PSFILE=%TEMP%\comet_setup_%RANDOM%.ps1"
> "%PSFILE%" echo $content = Get-Content -Raw -LiteralPath '%TEMPLATE%'
>> "%PSFILE%" echo $content = $content -replace 'name_of_project: "my_project"', 'name_of_project: "%PROJECT%"'
>> "%PSFILE%" echo $content = $content -replace 'HALO_data_file: "your_export.csv"', 'HALO_data_file: "%CSV_NAME%"'
>> "%PSFILE%" echo Set-Content -LiteralPath '%RUNFOLDER%\config.yaml' -Value $content

powershell -NoProfile -ExecutionPolicy Bypass -File "%PSFILE%"
if errorlevel 1 (
    echo.
    echo ERROR: Failed to write config.yaml. You can still copy
    echo config.yaml.template into "%RUNFOLDER%" by hand and edit it.
    echo.
    del "%PSFILE%" >nul 2>&1
    pause
    exit /b 1
)
del "%PSFILE%" >nul 2>&1

echo.
echo ============================================================
echo   Done! Your run folder is ready:
echo   %RUNFOLDER%
echo ============================================================
echo.
echo Next steps:
echo   1. Open that folder's config.yaml in Notepad and double-check
echo      it - project name, PCA dims, clustering resolution.
echo   2. If you already have a metadata_markers.csv for this marker
echo      panel, copy it into that same folder now.
echo      If you don't, that's fine - the pipeline will generate one
echo      automatically the first time it runs, then pause so you can
echo      review it. Just run RUN_COMET_PIPELINE.BAT a second time
echo      once you've checked it.
echo   3. Run RUN_COMET_PIPELINE.BAT and paste this path in when asked:
echo      %RUNFOLDER%
echo.
pause
