@echo off
setlocal
title COMET Pipeline Launcher

rem ─── Move to the folder containing this .bat file (the pipeline repo root) ─
cd /d "%~dp0"

rem ─── Check Python is available ────────────────────────────────────────────
where python >nul 2>&1
if %errorlevel% neq 0 (
    echo.
    echo  ERROR: Python was not found on this computer.
    echo  Please install Python 3.x from https://www.python.org/
    echo  and make sure "Add Python to PATH" is checked during installation.
    echo.
    pause
    exit /b 1
)

rem ─── Launch the GUI ────────────────────────────────────────────────────────
python launch_pipeline_gui.py
