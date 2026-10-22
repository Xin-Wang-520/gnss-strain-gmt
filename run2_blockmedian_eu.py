#!/usr/bin/env python3
"""Run 2: blockmedian of Eurasia-fixed GNSS, without velocity interpolation.

python3 run2_blockmedian_eu.py
python3 run2_blockmedian_eu.py --spacing 15m
Default: 30 arcmin bins over 67/113/22/45, including the input buffer.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
from run1_plot_gnss import ROOT, read_gnss
from GNSS_analysis.gnss_processing import blockmedian_input, locate_gmt, gmt_environment, run_command
import numpy as np

DEFAULT_REGION = [67, 113, 22, 45]
DEFAULT_OUTPUT = ROOT / 'GNSS_analysis/output/run2_blockmedian_eu/bm30m'


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spacing', '--block-spacing', default='30m', help='GMT blockmedian -I; m means arcminutes')
    parser.add_argument('--region', nargs=4, type=float, default=DEFAULT_REGION)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--gmt')
    args = parser.parse_args()
    parts = args.spacing.split('/')
    if len(parts) > 2 or any(not re.fullmatch(r'(?:\d+(?:\.\d*)?|\.\d+)[dmsk]?', p) or float(p.rstrip('dmsk')) <= 0 for p in parts):
        parser.error('Use positive GMT block increments, e.g. 30m or 10k/10k')
    w, e, s, n = args.region
    if not (67 <= w < e <= 113 and 22 <= s < n <= 45):
        parser.error('Region must lie inside 67/113/22/45')
    label = 'bm' + args.spacing.replace('.', 'p').replace('/', 'x')
    if args.region != DEFAULT_REGION:
        label += '_' + hashlib.sha256(json.dumps(args.region).encode()).hexdigest()[:8]
    args.output_dir = (args.output_dir or ROOT / 'GNSS_analysis/output/run2_blockmedian_eu' / label).resolve()
    return args


def main():
    args = arguments()
    data, _ = read_gnss(ROOT / 'GNSS_analysis/data/raw', args.region, frame='eu')
    gmt, env = locate_gmt(args.gmt), gmt_environment()
    reduced, metadata = blockmedian_input(gmt, data, args.region, args.spacing, args.output_dir, env)
    destination = args.output_dir / 'gnss_blockmedian_eu.txt'
    np.savetxt(destination, reduced, fmt='%.12g', header='lon lat VE VN sigmaE sigmaN; mm/yr')
    def study_count(array):
        return int(((array[:, 0] >= 70) & (array[:, 0] <= 110) &
                    (array[:, 1] >= 25) & (array[:, 1] <= 42)).sum())
    record = {'script': Path(__file__).name, 'stage': 2, 'method': 'GMT component-wise unweighted blockmedian',
              'gmt_version': run_command([gmt, '--version'], env, capture_output=True, text=True).stdout.strip(),
              'frame': 'ITRF2014 Eurasia-fixed', 'input_region': args.region, 'blockmedian': metadata,
              'output_file': str(destination), 'output_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
              'original_study_count': study_count(data), 'block_study_count': study_count(reduced),
              'columns': ['lon', 'lat', 'VE', 'VN', 'sigmaE', 'sigmaN'],
              'velocity_units': 'mm/yr', 'next_step': 'run3_gpsgridder_eu.py --input-dir ' + str(args.output_dir)}
    (args.output_dir / 'run2_summary.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'Run 2 blockmedian {args.spacing}: {len(data)} stations -> {len(reduced)} representative points.')
    print(f'Within 70/110/25/42: {study_count(data)} stations -> {study_count(reduced)} representatives.')
    print(f'Output: {destination}')


if __name__ == '__main__':
    main()
