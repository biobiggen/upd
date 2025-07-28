#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UPD Results Visualization
-------------------------
This script provides additional visualizations for UPD detection results.
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch


def load_results(results_dir):
    """
    Load UPD detection results from the output directory.
    """
    summary_file = os.path.join(results_dir, "upd_summary.tsv")
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    
    summary_df = pd.read_csv(summary_file, sep='\t')
    print(f"Loaded summary data for {len(summary_df)} chromosomes")
    
    return summary_df


def create_ideogram(output_dir, summary_df):
    """
    Create a chromosome ideogram showing UPD regions.
    """
    # Check if we have any data
    if len(summary_df) == 0:
        print("No data available for ideogram visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No UPD data available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "upd_ideogram.png"), dpi=150)
        plt.close()
        return os.path.join(output_dir, "upd_ideogram.png")
    
    # Define chromosome lengths (approximate)
    chr_lengths = {
        '1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276,
        '5': 180915260, '6': 171115067, '7': 159138663, '8': 146364022,
        '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895,
        '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753,
        '17': 81195210, '18': 78077248, '19': 59128983, '20': 63025520,
        '21': 48129895, '22': 51304566, 'X': 155270560
    }
    
    # Define UPD type colors
    upd_colors = {
        'No UPD detected': 'green',
        'Partial Maternal UPD': 'orange',
        'Partial Paternal UPD': 'purple',
        'Complete Maternal Isodisomy': 'darkred',
        'Complete Maternal Heterodisomy': 'red',
        'Mixed Maternal Iso/Heterodisomy': 'salmon',
        'Complete Paternal Isodisomy': 'darkblue',
        'Complete Paternal Heterodisomy': 'blue',
        'Mixed Paternal Iso/Heterodisomy': 'lightblue',
        'Insufficient data': 'gray'
    }
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Sort chromosomes
    chromosomes = sorted(summary_df['Chromosome'].unique(), 
                         key=lambda x: (0 if x.isdigit() else 1, int(x) if x.isdigit() else ord(x)))
    
    # Plot chromosomes as bars
    y_positions = []
    y_labels = []
    
    for i, chrom in enumerate(chromosomes):
        y_pos = i * 2
        y_positions.append(y_pos)
        y_labels.append(chrom)
        
        # Get UPD type for this chromosome
        chrom_data = summary_df[summary_df['Chromosome'] == chrom]
        if len(chrom_data) > 0:
            upd_type = chrom_data.iloc[0]['UPD_Type']
            color = upd_colors.get(upd_type, 'gray')
        else:
            color = 'gray'
        
        # Plot chromosome
        length = chr_lengths.get(chrom, 100000000)  # Default length if not in dictionary
        plt.barh(y_pos, length, height=1.0, color=color, alpha=0.7)
        
        # Add chromosome label
        plt.text(-10000000, y_pos, chrom, ha='right', va='center', fontsize=10)
    
    # Create legend
    legend_elements = [Patch(facecolor=color, alpha=0.7, label=upd_type) 
                      for upd_type, color in upd_colors.items() 
                      if upd_type in summary_df['UPD_Type'].values]
    
    plt.legend(handles=legend_elements, loc='upper right', fontsize=8)
    
    # Set plot properties
    plt.yticks([])  # Hide y-axis ticks
    plt.xlabel('Genomic Position')
    plt.title('Chromosome UPD Status')
    plt.xlim(-20000000, 250000000)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(results_dir, "upd_ideogram.png"), dpi=150)
    plt.close()


def create_heatmap(summary_df, output_dir):
    """
    Create a heatmap showing UPD scores and heterozygosity across chromosomes.
    """
    # Check if we have any data
    if len(summary_df) == 0:
        print("No data available for heatmap visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No UPD data available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "upd_heatmap.png"), dpi=150)
        plt.close()
        return os.path.join(output_dir, "upd_heatmap.png")
    
    # Check for NaN values in UPD_Score and Heterozygosity
    valid_data = summary_df.dropna(subset=['UPD_Score', 'Heterozygosity'])
    if len(valid_data) == 0:
        print("No valid numeric data available for heatmap visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No valid UPD metrics available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        plt.savefig(os.path.join(output_dir, "upd_heatmap.png"), dpi=150)
        plt.close()
        return os.path.join(output_dir, "upd_heatmap.png")
    
    # Prepare data
    heatmap_data = valid_data.copy()
    
    # Create UPD type numeric mapping for color gradient
    upd_type_map = {
        'No UPD detected': 0,
        'Partial Maternal UPD': 1,
        'Partial Paternal UPD': 2,
        'Complete Maternal Isodisomy': 3,
        'Complete Maternal Heterodisomy': 4,
        'Mixed Maternal Iso/Heterodisomy': 5,
        'Complete Paternal Isodisomy': 6,
        'Complete Paternal Heterodisomy': 7,
        'Mixed Paternal Iso/Heterodisomy': 8,
        'Insufficient data': -1
    }
    
    heatmap_data['UPD_Type_Numeric'] = heatmap_data['UPD_Type'].map(upd_type_map)
    
    # Sort chromosomes
    heatmap_data['Chromosome_Numeric'] = heatmap_data['Chromosome'].apply(
        lambda x: int(x) if x.isdigit() else (23 if x == 'X' else 24)
    )
    heatmap_data = heatmap_data.sort_values('Chromosome_Numeric')
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Create heatmap data
    df_heatmap = heatmap_data.pivot_table(
        index='Chromosome', 
        values=['UPD_Score', 'Heterozygosity', 'UPD_Type_Numeric'],
        aggfunc='first'
    )
    
    # Define colormap
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    
    # Plot heatmap
    sns.heatmap(df_heatmap, cmap=cmap, annot=True, fmt=".3f",
                linewidths=0.5, center=0)
    
    plt.title('UPD Detection Metrics by Chromosome')
    plt.tight_layout()
    
    # Save figure
    heatmap_file = os.path.join(output_dir, "upd_heatmap.png")
    plt.savefig(heatmap_file, dpi=150)
    plt.close()
    
    return heatmap_file


def create_detailed_report(results_dir, summary_df):
    """
    Create a detailed HTML report of UPD detection results.
    """
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>UPD Detection Results</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            h1 { color: #333366; }
            table { border-collapse: collapse; width: 100%; margin-top: 20px; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            tr:nth-child(even) { background-color: #f9f9f9; }
            .upd-normal { background-color: #c8e6c9; }
            .upd-maternal { background-color: #ffcdd2; }
            .upd-paternal { background-color: #bbdefb; }
            .upd-insufficient { background-color: #eeeeee; }
            .summary-section { margin-top: 30px; }
            .image-container { margin-top: 20px; text-align: center; }
            .image-container img { max-width: 100%; height: auto; margin-bottom: 10px; }
        </style>
    </head>
    <body>
        <h1>UPD Detection Results</h1>
        
        <div class="summary-section">
            <h2>Summary</h2>
            <p>Total chromosomes analyzed: {total_chromosomes}</p>
            <p>Chromosomes with UPD: {upd_chromosomes}</p>
            <ul>
                <li>Maternal UPD: {maternal_upd}</li>
                <li>Paternal UPD: {paternal_upd}</li>
            </ul>
        </div>
        
        <div class="summary-section">
            <h2>Chromosome Results</h2>
            <table>
                <tr>
                    <th>Chromosome</th>
                    <th>UPD Type</th>
                    <th>Heterozygosity</th>
                    <th>UPD Score</th>
                    <th>SNP Count</th>
                </tr>
                {table_rows}
            </table>
        </div>
        
        <div class="image-container">
            <h2>Visualizations</h2>
            <div>
                <img src="upd_summary.png" alt="UPD Summary">
                <p>UPD Summary Visualization</p>
            </div>
            <div>
                <img src="upd_ideogram.png" alt="UPD Ideogram">
                <p>Chromosome UPD Status</p>
            </div>
            <div>
                <img src="upd_heatmap.png" alt="UPD Heatmap">
                <p>UPD Detection Metrics Heatmap</p>
            </div>
        </div>
    </body>
    </html>
    """
    
    # Count UPD chromosomes
    total_chromosomes = len(summary_df)
    upd_chromosomes = len(summary_df[summary_df['UPD_Type'] != 'No UPD detected'])
    maternal_upd = len(summary_df[summary_df['UPD_Type'].str.contains('Maternal')])
    paternal_upd = len(summary_df[summary_df['UPD_Type'].str.contains('Paternal')])
    
    # Generate table rows
    table_rows = ""
    for _, row in summary_df.iterrows():
        upd_class = "upd-normal"
        if "Maternal" in row['UPD_Type']:
            upd_class = "upd-maternal"
        elif "Paternal" in row['UPD_Type']:
            upd_class = "upd-paternal"
        elif "Insufficient" in row['UPD_Type']:
            upd_class = "upd-insufficient"
            
        table_rows += f"""
        <tr class="{upd_class}">
            <td>{row['Chromosome']}</td>
            <td>{row['UPD_Type']}</td>
            <td>{row['Heterozygosity']:.3f}</td>
            <td>{row['UPD_Score']:.3f}</td>
            <td>{row['SNP_Count']}</td>
        </tr>
        """
    
    # Fill in the template
    html_content = html_content.format(
        total_chromosomes=total_chromosomes,
        upd_chromosomes=upd_chromosomes,
        maternal_upd=maternal_upd,
        paternal_upd=paternal_upd,
        table_rows=table_rows
    )
    
    # Write HTML file
    with open(os.path.join(results_dir, "upd_report.html"), 'w') as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description='Visualize UPD detection results')
    parser.add_argument('--results', '-r', default='upd_results', help='Results directory')
    
    args = parser.parse_args()
    
    # Load results
    summary_df = load_results(args.results)
    
    # Create visualizations
    create_ideogram(args.results, summary_df)
    heatmap_file = create_heatmap(summary_df, args.results)
    create_detailed_report(args.results, summary_df)
    
    print(f"Visualizations created in {args.results}")
    print(f"Detailed report: {os.path.join(args.results, 'upd_report.html')}")


if __name__ == "__main__":
    main()