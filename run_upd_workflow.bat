@echo off
echo UPD Detection Workflow Runner
echo ===========================

REM Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Error: Python is not installed or not in PATH
    echo Please install Python 3.6+ and try again
    pause
    exit /b 1
)

REM Check for required packages
echo Checking required packages...
python -c "import numpy, pandas, matplotlib, scipy, sklearn" >nul 2>&1
if %errorlevel% neq 0 (
    echo Installing required packages...
    pip install numpy pandas matplotlib scipy scikit-learn
)

echo.
echo Choose an option:
echo 1. Generate sample data and run UPD detection
echo 2. Use existing BAF data file
echo.

set /p option="Enter option (1 or 2): "

if "%option%"=="1" (
    echo.
    echo Generating sample data with UPD...
    set /p complete="Enter chromosomes with complete UPD (comma-separated, default: 7): " || set complete=7
    set /p partial="Enter chromosomes with partial UPD (comma-separated, default: 15): " || set partial=15
    set /p fetal_fraction="Enter estimated fetal fraction (default: 0.1): " || set fetal_fraction=0.1
    
    python run_upd_workflow.py --generate-sample --complete-upd %complete% --partial-upd %partial% --fetal-fraction %fetal_fraction%
) else if "%option%"=="2" (
    echo.
    set /p input_file="Enter path to BAF data file: "
    set /p fetal_fraction="Enter estimated fetal fraction (default: 0.1): " || set fetal_fraction=0.1
    
    if not exist "%input_file%" (
        echo Error: File not found: %input_file%
        pause
        exit /b 1
    )
    
    python run_upd_workflow.py --input "%input_file%" --fetal-fraction %fetal_fraction%
) else (
    echo Invalid option selected
    pause
    exit /b 1
)

echo.
echo Analysis complete!
echo Results are saved in the upd_results directory
echo.

pause