#!/bin/bash

echo "UPD Detection Workflow Runner"
echo "==========================="

# Check if Python is installed
if ! command -v python3.9 &> /dev/null; then
    echo "Error: Python 3.9 is not installed or not in PATH"
    echo "Please install Python 3.9+ and try again"
    exit 1
fi

# Check for required packages
echo "Checking required packages..."
if ! python3.9 -c "import numpy, pandas, matplotlib, scipy, sklearn" &> /dev/null; then
    echo "Installing required packages..."
    pip3.9 install -r requirements.txt
fi

echo ""
echo "Choose an option:"
echo "1. Generate sample data and run UPD detection"
echo "2. Use existing BAF data file"
echo ""

read -p "Enter option (1 or 2): " option

if [ "$option" == "1" ]; then
    echo ""
    echo "Generating sample data with UPD..."
    read -p "Enter chromosomes with complete UPD (comma-separated, default: 7): " complete
    complete=${complete:-7}
    
    read -p "Enter chromosomes with partial UPD (comma-separated, default: 15): " partial
    partial=${partial:-15}
    
    read -p "Enter estimated fetal fraction (default: 0.1): " fetal_fraction
    fetal_fraction=${fetal_fraction:-0.1}
    
    python3.9 run_upd_workflow.py --generate-sample --complete-upd "$complete" --partial-upd "$partial" --fetal-fraction "$fetal_fraction"
elif [ "$option" == "2" ]; then
    echo ""
    read -p "Enter path to BAF data file: " input_file
    read -p "Enter estimated fetal fraction (default: 0.1): " fetal_fraction
    fetal_fraction=${fetal_fraction:-0.1}
    
    if [ ! -f "$input_file" ]; then
        echo "Error: File not found: $input_file"
        exit 1
    fi
    
    python3.9 run_upd_workflow.py --input "$input_file" --fetal-fraction "$fetal_fraction"
else
    echo "Invalid option selected"
    exit 1
fi

echo ""
echo "Analysis complete!"
echo "Results are saved in the upd_results directory"
echo ""

# Run visualization
echo "Creating additional visualizations..."
python3.9 visualize_results.py

echo "Done!"