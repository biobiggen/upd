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
Command-line entry point for the UPD calculation tool

Subcommands:
    single - single-sample UPD calculation
    batch  - parallel batch calculation
    report - aggregate JSON results into a CSV report

Usage:
    python -m upd_tool.cli single -i sample.readslist -o result.json -p probe.bed
    python -m upd_tool.cli batch -i ./readslist_dir/ -o ./results/ --threads 8 -p probe.bed
    python -m upd_tool.cli report --results ./results/ -o upd_report.csv
"""

import argparse
import csv
import json
import os
import re
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .core import PROBE_VERSIONS, UPDCalculator, convert_numpy_types
from .regions import REPORT_REGIONS

# Default readlist file suffix
READLIST_SUFFIX = '_consensus.mapped.clipped.snp.reads.list'

# Reserved JSON key for non-region information (underscore prefix avoids
# clashing with region names)
FF_INFO_KEY = '_ff_info'


# ============================================================================
# Shared helpers
# ============================================================================
def extract_sample_name(filename: str) -> str:
    """Extract the sample name from a readlist or JSON file name."""
    basename = os.path.basename(filename)
    basename = re.sub(r'_upd_results\.json$', '', basename)
    basename = re.sub(
        r'_consensus\.mapped\.clipped\.snp\.reads\.list$', '', basename
    )
    basename = re.sub(r'\.snp\.reads\.list$', '', basename)
    return basename


def _ensure_dir(path: str) -> None:
    """Create the directory if it does not exist."""
    if path and not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


# ============================================================================
# single subcommand
# ============================================================================
def run_single(args: argparse.Namespace) -> int:
    """Single-sample UPD calculation."""
    if not os.path.exists(args.input):
        print(f'Error: input file {args.input} does not exist', file=sys.stderr)
        return 1

    if args.probe_file and not os.path.exists(args.probe_file):
        print(f'Error: probe file {args.probe_file} does not exist',
              file=sys.stderr)
        return 1

    try:
        calculator = UPDCalculator(
            probe_version=args.probe_version,
            probe_file=args.probe_file,
            test_code=args.test_code,
        )
        calculator.load_readlist(args.input)

        print('Calculating UPD...')
        upd_results = calculator.predict_upd_hmm()
        upd_results = convert_numpy_types(upd_results)
        upd_results[FF_INFO_KEY] = convert_numpy_types(calculator.ff_info)

        _ensure_dir(os.path.dirname(args.output))
        with open(args.output, 'w', encoding='utf-8') as f:
            json.dump(upd_results, f, indent=2, ensure_ascii=False)

        print(f'UPD calculation finished, results saved to {args.output}')

        ff_info = calculator.ff_info
        print(f"\nFetal fraction: {ff_info['ff_used']:.4f} "
              f"(method: {ff_info['ff_method']})")
        if ff_info['ff_heterozygous'] is not None:
            print(f"  Maternal homozygous estimate: "
                  f"{ff_info['ff_homozygous']:.4f}")
            print(f"  Maternal heterozygous estimate: "
                  f"{ff_info['ff_heterozygous']:.4f}")
            print(f"  Ratio: {ff_info['ff_ratio']:.2f}")
            if ff_info['ff_method'] == 'heterozygous':
                print('  Note: ratio is clearly high, '
                      'genome-wide homozygosity suspected')
        else:
            print('  Too few maternal heterozygous sites for a cross-check')

        print('\nUPD result summary:')
        for region, result in upd_results.items():
            if region == FF_INFO_KEY:
                continue
            if result['status'] == 'success':
                print(f"  {region}: {result['final_state']} ({result['final_ratio']:.2%})")
            else:
                print(f"  {region}: {result['status']}")

        # Optional plotting
        if args.plot:
            from .plotting import plot_from_readlist
            _ensure_dir(os.path.dirname(args.plot))
            n_points = plot_from_readlist(
                calculator.readlist_df,
                args.plot,
                title=calculator.cli,
            )
            if n_points > 0:
                print(f'Scatter plot saved to: {args.plot} ({n_points} points)')
            else:
                print('Warning: no SNP sites available to plot', file=sys.stderr)

        return 0

    except Exception as e:
        print(f'Error while calculating UPD: {e}', file=sys.stderr)
        return 1


# ============================================================================
# batch subcommand
# ============================================================================
def _process_one_sample(
    task: Tuple[str, str, Optional[str], str, int, Optional[str]]
) -> Tuple[str, bool, str]:
    """Process a single sample (called from the process pool).

    Args:
        task: (readlist_path, json_output, probe_file, probe_version,
               test_code, plot_output)

    Returns:
        (sample_name, success, message)
    """
    readlist_path, json_output, probe_file, probe_version, test_code, plot_output = task
    sample_name = extract_sample_name(readlist_path)

    try:
        calculator = UPDCalculator(
            probe_version=probe_version,
            probe_file=probe_file,
            test_code=test_code,
        )
        calculator.load_readlist(readlist_path)
        upd_results = calculator.predict_upd_hmm()
        upd_results = convert_numpy_types(upd_results)
        upd_results[FF_INFO_KEY] = convert_numpy_types(calculator.ff_info)

        _ensure_dir(os.path.dirname(json_output))
        with open(json_output, 'w', encoding='utf-8') as f:
            json.dump(upd_results, f, indent=2, ensure_ascii=False)

        if plot_output:
            from .plotting import plot_from_readlist
            _ensure_dir(os.path.dirname(plot_output))
            plot_from_readlist(
                calculator.readlist_df, plot_output, title=sample_name
            )

        return sample_name, True, 'success'

    except Exception as e:
        return sample_name, False, str(e)


def _collect_readlist_files(input_dir: str, recursive: bool) -> List[Path]:
    """Collect readlist files."""
    base = Path(input_dir)
    pattern = '*.snp.reads.list'
    if recursive:
        return sorted(base.rglob(pattern))
    return sorted(base.glob(pattern))


def run_batch(args: argparse.Namespace) -> int:
    """Parallel batch UPD calculation."""
    if not os.path.isdir(args.input):
        print(f'Error: input directory {args.input} does not exist',
              file=sys.stderr)
        return 1

    if args.probe_file and not os.path.exists(args.probe_file):
        print(f'Error: probe file {args.probe_file} does not exist',
              file=sys.stderr)
        return 1

    files = _collect_readlist_files(args.input, not args.no_recursive)
    if not files:
        print(f'Error: no *.snp.reads.list files found in {args.input}',
              file=sys.stderr)
        return 1

    print(f'Found {len(files)} readlist files')

    _ensure_dir(args.output)
    plot_dir = args.plot_dir
    if plot_dir:
        _ensure_dir(plot_dir)

    tasks = []
    for fpath in files:
        sample_name = extract_sample_name(str(fpath))
        json_output = os.path.join(args.output, f'{sample_name}_upd_results.json')
        plot_output = (
            os.path.join(plot_dir, f'{sample_name}_upd_regions_pa_ratio.png')
            if plot_dir else None
        )
        tasks.append((
            str(fpath), json_output, args.probe_file,
            args.probe_version, args.test_code, plot_output,
        ))

    threads = max(1, args.threads)
    print(f'Processing in parallel with {threads} processes...\n')

    success_count = 0
    fail_count = 0

    with ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(_process_one_sample, task): task[0]
            for task in tasks
        }
        for future in as_completed(futures):
            sample_name, ok, message = future.result()
            if ok:
                success_count += 1
                if args.verbose:
                    print(f'  [OK] {sample_name}: {message}')
            else:
                fail_count += 1
                print(f'  [FAIL] {sample_name}: {message}', file=sys.stderr)

    print(f"\n{'=' * 50}")
    print(f'Done: {success_count} succeeded, {fail_count} failed')
    print(f'Results directory: {args.output}')
    if plot_dir:
        print(f'Image directory: {plot_dir}')
    print(f"{'=' * 50}")

    return 0 if fail_count == 0 else 1


# ============================================================================
# report subcommand
# ============================================================================
def collect_upd_results(json_files: List[str]) -> Dict[str, Dict]:
    """Collect UPD results from multiple JSON files."""
    upd_results = {}
    for json_file in json_files:
        if not os.path.exists(json_file):
            continue
        try:
            sample_name = extract_sample_name(json_file)
            with open(json_file, 'r', encoding='utf-8') as f:
                upd_results[sample_name] = json.load(f)
        except Exception as e:
            print(f'Error while processing file {json_file}: {e}',
                  file=sys.stderr)
    return upd_results


def generate_report(upd_results: Dict[str, Dict], output_file: str) -> int:
    """Generate the UPD CSV report.

    Returns:
        Number of data rows written
    """
    header_parts = ['Sample_Name']
    for region in REPORT_REGIONS:
        header_parts.append(f'{region}_UPD')
        header_parts.append(f'{region}_Ratio')
    header_parts += ['FF_Used', 'FF_Method', 'FF_Ratio']

    _ensure_dir(os.path.dirname(output_file))

    row_count = 0
    with open(output_file, 'w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header_parts)

        for sample_name in sorted(upd_results.keys()):
            data = upd_results[sample_name]
            result_parts = [sample_name]

            for region in REPORT_REGIONS:
                if region in data:
                    result = data[region]
                    status = result.get('status', 'unknown')
                    final_state = (
                        result.get('final_state', 'N/A')
                        if status == 'success' else status
                    )
                    final_ratio = (
                        result.get('final_ratio', 'N/A')
                        if status == 'success' else 'N/A'
                    )

                    result_parts.append(final_state)
                    if isinstance(final_ratio, (float, int)):
                        result_parts.append(f'{final_ratio:.4f}')
                    else:
                        result_parts.append(str(final_ratio))
                else:
                    result_parts.append('N/A')
                    result_parts.append('N/A')

            ff_info = data.get(FF_INFO_KEY) or {}
            for key, fmt in (('ff_used', '.4f'), ('ff_method', None),
                             ('ff_ratio', '.2f')):
                value = ff_info.get(key)
                if value is None:
                    result_parts.append('N/A')
                elif fmt and isinstance(value, (float, int)):
                    result_parts.append(format(value, fmt))
                else:
                    result_parts.append(str(value))

            writer.writerow(result_parts)
            row_count += 1

    return row_count


def run_report(args: argparse.Namespace) -> int:
    """Aggregate JSON results into a CSV report."""
    json_files: List[str] = []

    if args.results:
        if not os.path.isdir(args.results):
            print(f'Error: results directory {args.results} does not exist',
                  file=sys.stderr)
            return 1
        json_files = [
            str(p) for p in sorted(Path(args.results).rglob('*_upd_results.json'))
        ]
    elif args.upd_jsons:
        json_files = args.upd_jsons
    else:
        print('Error: either --results directory or --upd-jsons file list '
              'must be given', file=sys.stderr)
        return 1

    if not json_files:
        print('Error: no UPD result JSON files found', file=sys.stderr)
        return 1

    print(f'Found {len(json_files)} UPD result files')

    upd_results = collect_upd_results(json_files)
    if not upd_results:
        print('Error: no valid results could be loaded', file=sys.stderr)
        return 1

    row_count = generate_report(upd_results, args.output)
    print(f'Report generated: {args.output} ({row_count} samples)')
    return 0


# ============================================================================
# Argument parsing
# ============================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='upd_tool',
        description='UPD (Uniparental Disomy) calculation tool',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    probe_versions = list(PROBE_VERSIONS.keys())

    # --- single ---
    p_single = subparsers.add_parser('single',
                                     help='single-sample UPD calculation')
    p_single.add_argument('-i', '--input', required=True,
                          help='path to the input snp.readslist.txt file')
    p_single.add_argument('-o', '--output', required=True,
                          help='output result file path (JSON)')
    p_single.add_argument('-p', '--probe-file', default=None,
                          help='path to the probe file')
    p_single.add_argument('--probe-version', default='NIPT3V4',
                          choices=probe_versions,
                          help='probe version (default: NIPT3V4)')
    p_single.add_argument('-t', '--test-code', type=int, default=1,
                          help='test code type (default: 1)')
    p_single.add_argument('--plot', default=None,
                          help='also generate a scatter plot at this PNG path')
    p_single.set_defaults(func=run_single)

    # --- batch ---
    p_batch = subparsers.add_parser('batch',
                                    help='parallel batch UPD calculation')
    p_batch.add_argument('-i', '--input', required=True,
                         help='directory containing *.snp.reads.list files')
    p_batch.add_argument('-o', '--output', required=True,
                         help='output directory for JSON results')
    p_batch.add_argument('-p', '--probe-file', default=None,
                         help='path to the probe file')
    p_batch.add_argument('--probe-version', default='NIPT3V4',
                         choices=probe_versions,
                         help='probe version (default: NIPT3V4)')
    p_batch.add_argument('-t', '--test-code', type=int, default=1,
                         help='test code type (default: 1)')
    p_batch.add_argument('--threads', type=int, default=4,
                         help='number of parallel processes (default: 4)')
    p_batch.add_argument('--plot-dir', default=None,
                         help='also generate scatter plots in this directory')
    p_batch.add_argument('--no-recursive', action='store_true',
                         help='do not scan subdirectories (default: recursive)')
    p_batch.add_argument('-v', '--verbose', action='store_true',
                         help='print the result of every sample')
    p_batch.set_defaults(func=run_batch)

    # --- report ---
    p_report = subparsers.add_parser(
        'report', help='aggregate JSON results into a CSV report')
    p_report.add_argument('--results', default=None,
                          help='directory holding UPD result JSON files '
                               '(scanned recursively)')
    p_report.add_argument('--upd-jsons', nargs='+', default=None,
                          help='list of UPD result JSON files')
    p_report.add_argument('-o', '--output', required=True,
                          help='output CSV report path')
    p_report.set_defaults(func=run_report)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == '__main__':
    main()
