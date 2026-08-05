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
UPD (Uniparental Disomy) calculation tool

A standalone UPD analysis tool for NIPT SNP readlists, supporting both
single-sample and batch computation.

Modules:
    core     - core UPD computation (UPDCalculator)
    regions  - UPD target region coordinates
    plotting - pA_Ratio scatter plotting
    cli      - command-line entry point (single/batch/report)
"""

__version__ = '2.0.0'

from .core import UPDCalculator, PROBE_VERSIONS
from .regions import UPD_TARGETS_HMM, UPD_TARGETS_PLOT, REPORT_REGIONS

__all__ = [
    'UPDCalculator',
    'PROBE_VERSIONS',
    'UPD_TARGETS_HMM',
    'UPD_TARGETS_PLOT',
    'REPORT_REGIONS',
]
