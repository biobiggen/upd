#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright (C) 2024 biobiggen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
pA_Ratio scatter plotting for UPD regions

Accepts an already loaded DataFrame to avoid re-parsing the readlist file.
"""

import os

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from typing import Dict, Optional

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from .regions import PLOT_ORDER, UPD_TARGETS_PLOT

# Red / green / blue colour cycle
COLORS = ['#FF0000', '#00FF00', '#0000FF']


def filter_plot_regions(readlist_df: pd.DataFrame) -> pd.DataFrame:
    """Select SNP sites by plotting region and tag SNP_Tag with the region name.

    Args:
        readlist_df: the loaded readlist DataFrame

    Returns:
        Concatenated DataFrame with a SNP_Tag column holding the region name
    """
    filtered_dfs = []
    for region_name, region_info in UPD_TARGETS_PLOT.items():
        chr_filter = readlist_df['#Chr'] == region_info['chr']
        pos_filter = (
            (readlist_df['Pos'] >= region_info['start'])
            & (readlist_df['Pos'] <= region_info['end'])
        )
        region_df = readlist_df.loc[chr_filter & pos_filter, :].copy()
        if not region_df.empty:
            region_df['SNP_Tag'] = region_name
            filtered_dfs.append(region_df)

    if filtered_dfs:
        return pd.concat(filtered_dfs, ignore_index=True)
    return pd.DataFrame()


def plot_upd_regions(
    data: pd.DataFrame,
    output_file: str,
    title: Optional[str] = None,
) -> None:
    """Plot pA_Ratio scatter plot for UPD target regions.

    Args:
        data: DataFrame processed by :func:`filter_plot_regions`
        output_file: output image file path (PNG)
        title: plot title (derived from output filename by default)
    """
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    plt.figure(figsize=(15, 10))

    pos_slid = 0
    x_ticks = []
    x_labels = []

    for i, region_name in enumerate(PLOT_ORDER):
        tmp_df = data.loc[data.SNP_Tag == region_name]
        region_pa_ratio = list(tmp_df.pA_Ratio)

        x_positions = range(pos_slid, len(region_pa_ratio) + pos_slid)
        plt.scatter(
            x_positions,
            region_pa_ratio,
            label=f'{region_name} (n={len(region_pa_ratio)})',
            alpha=0.7,
            s=8,
            color=COLORS[i % 3],
        )

        if len(region_pa_ratio) > 0:
            x_ticks.append(pos_slid + len(region_pa_ratio) // 2)
            x_labels.append(region_name)

        pos_slid += len(region_pa_ratio)

    if title is None:
        basename = os.path.basename(output_file)
        title = '_'.join(basename.split('_')[:3])

    plt.title(f'UPD pA_Ratio Distribution {title}', fontsize=16)
    plt.ylabel('pA_Ratio')
    plt.ylim(0, 1)
    plt.xticks(x_ticks, x_labels, rotation=0, ha='center')
    plt.axhline(y=0.5, color='black', linestyle='--', alpha=0.5)

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()


def plot_from_readlist(
    readlist_df: pd.DataFrame,
    output_file: str,
    title: Optional[str] = None,
) -> int:
    """Plot directly from a readlist DataFrame.

    Args:
        readlist_df: the loaded readlist DataFrame
        output_file: output image path
        title: plot title

    Returns:
        Total number of plotted SNP sites
    """
    data = filter_plot_regions(readlist_df)
    if data.empty:
        return 0
    plot_upd_regions(data, output_file, title)
    return len(data)
