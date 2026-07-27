@echo off
setlocal enabledelayedexpansion
title COMET Pipeline Launcher

echo ============================================================
echo   COMET Spatial Pipeline
echo ============================================================
echo.
echo Your run folder (on the network drive) must already contain:
echo   - config.yaml           (your edited settings for this run)
echo   - your input CSV file   (path matching what's in config.yaml)
echo.
echo Don't have a run folder yet? Close this window and double-click
echo SETUP_NEW_RUN.BAT instead - it builds one for you.
echo.
echo Example: Z:\Yves\Run_2026_06_30
echo.

set "RUNFOLDER="
set /p RUNFOLDER=Enter the full path to your run folder:

rem Strip any surrounding quotes the user may have pasted in.
set RUNFOLDER=%RUNFOLDER:"=%

if "%RUNFOLDER%"=="" (
    echo.
    echo ERROR: No folder path entered.
    echo.
    pause
    exit /b 1
)

if not exist "%RUNFOLDER%\" (
    echo.
    echo ERROR: Could not find the folder:
    echo   "%RUNFOLDER%"
    echo Check the path and try again. If you copy-pasted it from File
    echo Explorer's address bar, make sure it points at the folder
    echo itself, not a file inside it.
    echo.
    pause
    exit /b 1
)

if not exist "%RUNFOLDER%\config.yaml" (
    echo.
    echo ERROR: No config.yaml found in "%RUNFOLDER%"
    echo Copy config.yaml.template from the pipeline_code folder into your
    echo run folder, edit it, save it as config.yaml, then run this again.
    echo ^(Or use SETUP_NEW_RUN.BAT to do this automatically.^)
    echo.
    pause
    exit /b 1
)

echo.
echo Checking Docker Desktop is running ...
docker info >nul 2>&1
if errorlevel 1 (
    echo.
    echo ERROR: Docker Desktop does not appear to be running.
    echo Start Docker Desktop from the Start Menu, wait until it says
    echo "Docker Desktop is running", then run this script again.
    echo.
    pause
    exit /b 1
)

set RUN_DIR=%RUNFOLDER%
cd /d "C:\Users\Yves\Comet_Data\Comet_Pipeline"

echo.
echo Starting pipeline for: %RUNFOLDER%
echo.
docker compose up --build

echo.
echo ============================================================
echo Pipeline finished or stopped.
echo If you saw an annotation step, did you finish it at http://localhost:8501 ?
echo Your results should be in: %RUNFOLDER%
echo ============================================================
pause
