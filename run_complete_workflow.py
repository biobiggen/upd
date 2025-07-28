#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Complete UPD Detection Workflow Runner
-------------------------------------
This script runs the complete workflow for detecting UPD in fetal DNA from maternal plasma:
1. Generate sample data or use provided BAF data
2. Run UPD detection with HMM
3. Generate visualizations and reports
"""

import os
import argparse
import subprocess
import time


def run_command(command, description=None):
    """Run a shell command and print output"""
    if description:
        print(f"\n=== {description} ===")
    
    print(f"Running: {' '.join(command)}")
    start_time = time.time()
    
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    
    # Print output in real-time
    for line in process.stdout:
        print(line.strip())
    
    # Wait for process to complete
    process.wait()
    
    # Print any errors
    if process.returncode != 0:
        print("Error output:")
        for line in process.stderr:
            print(line.strip())
        raise RuntimeError(f"Command failed with return code {process.returncode}")
    
    elapsed_time = time.time() - start_time
    print(f"Completed in {elapsed_time:.2f} seconds")
    
    return process.returncode


def main():
    parser = argparse.ArgumentParser(description='Run complete UPD detection workflow')
    parser.add_argument('--input', '-i', help='Input BAF file (if not provided, sample data will be generated)')
    parser.add_argument('--output', '-o', default='upd_results', help='Output directory')
    parser.add_argument('--generate-sample', '-g', action='store_true', help='Generate sample data')
    parser.add_argument('--complete-upd', default='7', help='Chromosomes with complete UPD (for sample data)')
    parser.add_argument('--partial-upd', default='15', help='Chromosomes with partial UPD (for sample data)')
    parser.add_argument('--fetal-fraction', '-f', type=float, default=0.1,
                        help='Estimated fetal DNA fraction in maternal plasma (default: 0.1)')
    parser.add_argument('--depth', '-d', type=int, default=100,
                        help='Estimated sequencing depth (default: 100)')
    parser.add_argument('--python', default='python3.9', help='Python executable to use')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Generate sample data if needed
    input_file = args.input
    if args.generate_sample or input_file is None:
        sample_file = os.path.join(args.output, 'sample_data.tsv')
        print(f"\n=== Generating sample data ===")
        print(f"Complete UPD chromosomes: {args.complete_upd}")
        print(f"Partial UPD chromosomes: {args.partial_upd}")
        
        generate_cmd = [
            args.python, 'generate_sample_data.py',
            '--output', sample_file,
            '--complete-upd', args.complete_upd,
            '--partial-upd', args.partial_upd,
            '--fetal-fraction', str(args.fetal_fraction)
        ]
        
        run_command(generate_cmd)
        input_file = sample_file
    
    # Run UPD detection
    upd_cmd = [
        args.python, 'upd_detection.py',
        '--input', input_file,
        '--output', args.output,
        '--fetal-fraction', str(args.fetal_fraction),
        '--depth', str(args.depth)
    ]
    
    run_command(upd_cmd, "Running UPD detection")
    
    # Generate additional visualizations
    vis_cmd = [
        args.python, 'visualize_upd_results.py',
        '--results', args.output
    ]
    
    run_command(vis_cmd, "Generating additional visualizations")
    
    print(f"\n=== Workflow completed successfully ===")
    print(f"Results saved to: {args.output}")
    print(f"Summary file: {os.path.join(args.output, 'upd_summary.tsv')}")
    print(f"Detailed report: {os.path.join(args.output, 'upd_report.html')}")


if __name__ == "__main__":
    main()