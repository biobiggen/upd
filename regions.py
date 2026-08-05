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
UPD target region coordinates

1. UPD_TARGETS_HMM  - narrow intervals used by the HMM computation
2. UPD_TARGETS_PLOT - wide intervals used for plotting
"""

from typing import Dict, List, Optional, Tuple


# ============================================================================
# HMM computation regions (narrow intervals)
# ============================================================================
# Format: region_name -> [chrom, start_pos, end_pos]
# 'chrRef' is the reference region, selected via SNP_Tag == 'chrRef' with no
# coordinate range
UPD_TARGETS_HMM: Dict[str, List] = {
    '6q24': ['chr6', 142448249, 145502506],
    '7q32': ['chr7', 129001172, 132163820],
    '11p15': ['chr11', 496924, 4384316],
    '14q32': ['chr14', 99229767, 103054425],
    '15q11q13': ['chr15', 23502875, 28268218],
    '20q13': ['chr20', 57354825, 60403476],
    'chrRef': ['chrRef', None, None],
}


# ============================================================================
# Plotting regions (wide intervals)
# ============================================================================
# Format: region_name -> {'chr': ..., 'start': ..., 'end': ...}
UPD_TARGETS_PLOT: Dict[str, Dict[str, object]] = {
    '6q24': {'chr': 'chr6', 'start': 138300000, 'end': 148500000},
    '7q32': {'chr': 'chr7', 'start': 127500000, 'end': 132900000},
    '11p15': {'chr': 'chr11', 'start': 0, 'end': 22000000},
    '14q32': {'chr': 'chr14', 'start': 89300000, 'end': 107043718},
    '15q11q13': {'chr': 'chr15', 'start': 19000000, 'end': 33400000},
    '20q13': {'chr': 'chr20', 'start': 43100000, 'end': 64444167},
}


# ============================================================================
# Plotting order (excluding chrRef)
# ============================================================================
PLOT_ORDER: Tuple[str, ...] = (
    '6q24',
    '7q32',
    '11p15',
    '14q32',
    '15q11q13',
    '20q13',
)


# ============================================================================
# Report region order (including chrRef)
# ============================================================================
REPORT_REGIONS: List[str] = [
    '6q24',
    '7q32',
    '11p15',
    '14q32',
    '15q11q13',
    '20q13',
    'chrRef',
]

