#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
UPD Results Visualization
------------------------
This script provides additional visualization options for UPD detection results.
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def load_summary(summary_file):
    """Load UPD summary results"""
    if not os.path.exists(summary_file):
        raise FileNotFoundError(f"Summary file not found: {summary_file}")
    
    return pd.read_csv(summary_file, sep='\t')


def create_heatmap(summary_df, output_dir):
    """Create a heatmap visualization of UPD results"""
    # Check if we have any data
    if len(summary_df) == 0:
        print("No data available for heatmap visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No UPD data available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        output_file = os.path.join(output_dir, "upd_heatmap.png")
        plt.savefig(output_file, dpi=150)
        plt.close()
        return output_file
    
    # Check for NaN values in UPD_Score and Heterozygosity
    valid_data = summary_df.dropna(subset=['UPD_Score', 'Heterozygosity'])
    if len(valid_data) == 0:
        print("No valid numeric data available for heatmap visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No valid UPD metrics available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        output_file = os.path.join(output_dir, "upd_heatmap.png")
        plt.savefig(output_file, dpi=150)
        plt.close()
        return output_file
    
    # Prepare data for heatmap
    chromosomes = valid_data['Chromosome'].tolist()
    
    # Create data matrix
    data = {
        'Heterozygosity': valid_data['Heterozygosity'].values,
        'UPD_Score': valid_data['UPD_Score'].values
    }
    
    # Create figure
    plt.figure(figsize=(12, 8))
    
    # Plot heatmap
    df_heatmap = pd.DataFrame(data, index=chromosomes)
    
    # Set up color mapping
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    
    # Plot heatmap
    sns.heatmap(df_heatmap, cmap=cmap, annot=True, fmt=".3f", 
                linewidths=.5, cbar_kws={"shrink": .8})
    
    plt.title("UPD Detection Results by Chromosome")
    plt.tight_layout()
    
    # Save figure
    output_file = os.path.join(output_dir, "upd_heatmap.png")
    plt.savefig(output_file, dpi=150)
    plt.close()
    
    return output_file


def create_circular_plot(summary_df, output_dir):
    """Create a circular plot of chromosomes with UPD highlighted"""
    # Check if we have any data
    if len(summary_df) == 0:
        print("No data available for circular plot visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No UPD data available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        output_file = os.path.join(output_dir, "upd_circular_plot.png")
        plt.savefig(output_file, dpi=150)
        plt.close()
        return output_file
    
    # Check for NaN values in Heterozygosity
    valid_data = summary_df.dropna(subset=['Heterozygosity'])
    if len(valid_data) == 0:
        print("No valid heterozygosity data available for circular plot visualization.")
        # Create a simple placeholder image
        plt.figure(figsize=(10, 6))
        plt.text(0.5, 0.5, "No valid heterozygosity data available for visualization", 
                 ha='center', va='center', fontsize=14)
        plt.axis('off')
        output_file = os.path.join(output_dir, "upd_circular_plot.png")
        plt.savefig(output_file, dpi=150)
        plt.close()
        return output_file
    
    try:
        import circlify
        import matplotlib.patches as patches
    except ImportError:
        print("circlify package not found. Installing...")
        import subprocess
        subprocess.check_call(["pip", "install", "circlify"])
        import circlify
        import matplotlib.patches as patches
    
    # Prepare data
    chromosomes = valid_data['Chromosome'].tolist()
    het_values = valid_data['Heterozygosity'].values
    upd_types = valid_data['UPD_Type'].values
    
    # Create size values (inverse of heterozygosity)
    sizes = 1 - het_values + 0.5  # Add 0.5 to ensure all chromosomes are visible
    
    # Create color values based on UPD type
    colors = []
    for upd_type in upd_types:
        if upd_type == "Complete UPD":
            colors.append("#FF5555")  # Red
        elif upd_type == "Partial UPD":
            colors.append("#FFAA55")  # Orange
        else:
            colors.append("#55AA55")  # Green
    
    # Create data for circlify
    data = [{"id": chr, "datum": size, "color": color} 
            for chr, size, color in zip(chromosomes, sizes, colors)]
    
    # Compute circle positions
    circles = circlify.circlify(
        data, 
        show_enclosure=False,
        target_enclosure=circlify.Circle(x=0, y=0, r=1)
    )
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Set equal aspect ratio
    ax.set_aspect('equal')
    ax.set_axis_off()
    
    # Plot circles
    for circle, label in zip(circles, data):
        x, y, r = circle.x, circle.y, circle.r
        ax.add_patch(patches.Circle((x, y), r, alpha=0.9, linewidth=2, 
                                    facecolor=label["color"]))
        plt.annotate(label["id"], (x, y), ha='center', va='center', 
                     fontsize=max(8, int(r*30)))
    
    # Add legend
    legend_elements = [
        patches.Patch(facecolor="#FF5555", edgecolor="w", label="Complete UPD"),
        patches.Patch(facecolor="#FFAA55", edgecolor="w", label="Partial UPD"),
        patches.Patch(facecolor="#55AA55", edgecolor="w", label="No UPD")
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=12)
    
    plt.title("Chromosome UPD Status", fontsize=16)
    
    # Save figure
    output_file = os.path.join(output_dir, "upd_circular_plot.png")
    plt.savefig(output_file, dpi=150)
    plt.close()
    
    return output_file


def create_interactive_html(summary_df, output_dir):
    """Create an interactive HTML report"""
    # Check if we have any data
    if len(summary_df) == 0:
        print("No data available for interactive HTML report.")
        # Create a simple HTML file
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>UPD Detection Results</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; text-align: center; }
                h1 { color: #333366; }
                .message { margin-top: 50px; font-size: 18px; color: #666; }
            </style>
        </head>
        <body>
            <h1>UPD Detection Results</h1>
            <div class="message">No UPD data available for visualization</div>
        </body>
        </html>
        """
        output_file = os.path.join(output_dir, "upd_interactive_report.html")
        with open(output_file, 'w') as f:
            f.write(html_content)
        return output_file
    
    # Check for NaN values
    valid_data = summary_df.dropna(subset=['UPD_Score', 'Heterozygosity'])
    if len(valid_data) == 0:
        print("No valid numeric data available for interactive HTML report.")
        # Create a simple HTML file
        html_content = """
        <!DOCTYPE html>
        <html>
        <head>
            <title>UPD Detection Results</title>
            <style>
                body { font-family: Arial, sans-serif; margin: 20px; text-align: center; }
                h1 { color: #333366; }
                .message { margin-top: 50px; font-size: 18px; color: #666; }
            </style>
        </head>
        <body>
            <h1>UPD Detection Results</h1>
            <div class="message">No valid UPD metrics available for visualization</div>
        </body>
        </html>
        """
        output_file = os.path.join(output_dir, "upd_interactive_report.html")
        with open(output_file, 'w') as f:
            f.write(html_content)
        return output_file
    
    try:
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError:
        print("plotly package not found. Installing...")
        import subprocess
        subprocess.check_call(["pip", "install", "plotly"])
        import plotly.express as px
        import plotly.graph_objects as go
    
    # Create bar chart
    fig = px.bar(
        valid_data, 
        x='Chromosome', 
        y='Heterozygosity',
        color='UPD_Type',
        hover_data=['UPD_Score', 'SNP_Count'],
        title='Chromosome Heterozygosity and UPD Status',
        color_discrete_map={
            'Complete UPD': '#FF5555',
            'Partial UPD': '#FFAA55',
            'No UPD detected': '#55AA55',
            'Complete Maternal Isodisomy': '#FF5555',
            'Complete Maternal Heterodisomy': '#FF7777',
            'Mixed Maternal Iso/Heterodisomy': '#FF9999',
            'Complete Paternal Isodisomy': '#5555FF',
            'Complete Paternal Heterodisomy': '#7777FF',
            'Mixed Paternal Iso/Heterodisomy': '#9999FF',
            'Partial Maternal UPD': '#FFAA55',
            'Partial Paternal UPD': '#55AAFF',
            'Insufficient data': '#AAAAAA'
        }
    )
    
    # Add threshold line
    fig.add_shape(
        type="line",
        x0=-0.5,
        y0=0.25,
        x1=len(summary_df)-0.5,
        y1=0.25,
        line=dict(color="red", width=2, dash="dash"),
    )
    
    # Add annotation for threshold
    fig.add_annotation(
        x=len(summary_df)-1,
        y=0.25,
        text="UPD Threshold",
        showarrow=False,
        yshift=10
    )
    
    # Update layout
    fig.update_layout(
        xaxis_title="Chromosome",
        yaxis_title="Heterozygosity Ratio",
        legend_title="UPD Status",
        font=dict(size=14)
    )
    
    # Save as HTML
    output_file = os.path.join(output_dir, "upd_interactive_report.html")
    fig.write_html(output_file)
    
    return output_file


def main():
    parser = argparse.ArgumentParser(description='Visualize UPD detection results')
    parser.add_argument('--input', '-i', default='upd_results/upd_summary.tsv', 
                        help='Input summary TSV file')
    parser.add_argument('--output', '-o', default=None, 
                        help='Output directory (defaults to same directory as input)')
    
    args = parser.parse_args()
    
    # Determine output directory
    if args.output:
        output_dir = args.output
    else:
        output_dir = os.path.dirname(args.input)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load summary data
    summary_df = load_summary(args.input)
    
    print(f"Loaded summary data with {len(summary_df)} chromosomes")
    
    # Create visualizations
    print("Creating heatmap visualization...")
    heatmap_file = create_heatmap(summary_df, output_dir)
    
    print("Creating circular plot visualization...")
    circular_file = create_circular_plot(summary_df, output_dir)
    
    print("Creating interactive HTML report...")
    html_file = create_interactive_html(summary_df, output_dir)
    
    print("\nVisualization complete!")
    print(f"Heatmap: {heatmap_file}")
    print(f"Circular plot: {circular_file}")
    print(f"Interactive report: {html_file}")


if __name__ == "__main__":
    main()