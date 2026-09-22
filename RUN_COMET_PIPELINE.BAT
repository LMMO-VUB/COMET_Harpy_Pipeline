@echo off
cd /d "%~dp0"

rem Use pythonw (no console window) if available, otherwise fall back to python.
where pythonw >nul 2>&1
if %errorlevel% equ 0 (
    start "" pythonw launch_pipeline_gui.py
) else (
    where python >nul 2>&1
    if %errorlevel% neq 0 (
        echo Python is not installed or not in PATH.
        echo Install from https://www.python.org/ with "Add to PATH" checked.
        pause
        exit /b 1
    )
    start "" python launch_pipeline_gui.py
)
