@echo off
setlocal
title COMET Pipeline Launcher

echo ============================================================
echo   COMET Spatial Pipeline
echo ============================================================
echo.
echo Your run folder (on the network drive) must already contain:
echo   - config.yaml           (your edited settings for this run)
echo   - your input CSV files  (paths matching what's in config.yaml)
echo.
echo Example: Z:\Yves\Run_2026_06_30
echo.

set /p RUNFOLDER=Enter the full path to your run folder:

if not exist "%RUNFOLDER%\config.yaml" (
    echo.
    echo ERROR: No config.yaml found in "%RUNFOLDER%"
    echo Copy config.yaml.template from the pipeline_code folder into your
    echo run folder, edit it, save it as config.yaml, then run this again.
    echo.
    pause
    exit /b 1
)

set RUN_DIR=%RUNFOLDER%
cd /d "C:\Users\Yves\Comet_Data\Comet_Pipeline"

docker compose up --build

echo.
echo ============================================================
echo Pipeline finished or stopped.
echo If you saw an annotation step, did you finish it at http://localhost:8501 ?
echo Your results should be in: %RUNFOLDER%
echo ============================================================
pause

