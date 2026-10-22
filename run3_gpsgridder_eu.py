#!/usr/bin/env python3
"""Run 3: gpsgridder on the blockmedian representatives produced by Run 2.

python3 run3_gpsgridder_eu.py -v
python3 run3_gpsgridder_eu.py --show-mode-counts
Defaults: run2 bm30m input, 25% modes, 20 arcmin output, Eurasia-fixed.
Identical completed fits are reused automatically; --rerun refits them.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import shlex
import time
from run1_plot_gnss import ROOT
from run2_blockmedian_eu import DEFAULT_OUTPUT as DEFAULT_INPUT
from GNSS_analysis.gnss_processing import (
    configuration, diagnostics, gmt_environment, increment_degrees, locate_gmt,
    report_used_modes, requested_mode_count, run_command, run_gpsgridder, token)
import numpy as np


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--region', type=float, nargs=4, default=[70, 110, 25, 42])
    parser.add_argument('--spacing', default='20m', help='gpsgridder output -I; independent of blockmedian -I')
    parser.add_argument('--retain-percent', type=float, default=25)
    parser.add_argument('--poisson', type=float)
    parser.add_argument('--fudge-km', type=float)
    parser.add_argument('--mask-distance-km', type=float, default=0)
    parser.add_argument('--gmt')
    parser.add_argument('--dpi', type=int, default=220)
    parser.add_argument('-v', '--verbose', action='store_true')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--rerun', action='store_true', help='Refit even if matching grids exist')
    mode.add_argument('--reuse-grids', action='store_true', help='Require an existing matching fit')
    parser.add_argument('--show-mode-counts', action='store_true', help='Only show counts; do not fit or draw')
    args = parser.parse_args()
    if not 0 < args.retain_percent <= 100 or (args.poisson is not None and not -1 <= args.poisson <= 1):
        parser.error('Use 0 < retained percentage <= 100 and -1 <= Poisson ratio <= 1')
    try:
        args.spacing_degree = increment_degrees(args.spacing)
    except (ValueError, IndexError):
        parser.error('Use angular output spacing, e.g. 20m')
    if not np.isfinite(args.spacing_degree) or args.spacing_degree <= 0 or args.dpi <= 0 or not np.isfinite(args.mask_distance_km) or args.mask_distance_km < 0:
        parser.error('Use positive spacing/dpi and nonnegative mask distance')
    if args.fudge_km is not None and (not np.isfinite(args.fudge_km) or args.fudge_km <= 0):
        parser.error('Distance offset must be positive and finite')
    w, e, s, n = args.region
    if not np.isfinite(args.region).all() or not (w < e and s < n) or any(not np.isclose(length / args.spacing_degree, round(length / args.spacing_degree)) for length in (e-w, n-s)):
        parser.error('Output region dimensions must be positive integer multiples of spacing')
    args.input_dir = args.input_dir.resolve()
    return args


def load_inputs(args):
    source = args.input_dir / 'run2_summary.json'
    if not source.is_file():
        raise FileNotFoundError(f'Run 2 input missing: {source}. Run run2_blockmedian_eu.py first.')
    manifest = json.loads(source.read_text())
    path = args.input_dir / 'gnss_blockmedian_eu.txt'
    if hashlib.sha256(path.read_bytes()).hexdigest() != manifest['output_sha256']:
        raise ValueError('Run 2 representative file changed since its summary was written')
    if manifest['frame'] != 'ITRF2014 Eurasia-fixed':
        raise ValueError('Expected Eurasia-fixed GNSS')
    data = np.loadtxt(path, ndmin=2)
    metadata = manifest['blockmedian']
    if data.shape != (metadata['occupied_block_count'], 6) or not np.isfinite(data).all() or np.any(data[:, 4:] <= 0):
        raise ValueError('Invalid six-column representative data')
    args.input_region = manifest['input_region']
    args.block_spacing = metadata['gmt_increment']
    iw, ie, iss, inn = args.input_region
    w, e, s, n = args.region
    if not (iw <= w < e <= ie and iss <= s < n <= inn):
        raise ValueError('Output region must lie inside the Run 2 input region')
    label = f'bm{args.block_spacing.replace(".", "p").replace("/", "x")}_c{token(args.retain_percent)}_i{args.spacing.replace(".", "p")}_pixel'
    if args.poisson is not None:
        label += '_nu' + token(args.poisson)
    if args.fudge_km is not None:
        label += '_fd' + token(args.fudge_km)
    if args.region != [70, 110, 25, 42] or args.input_region != [67, 113, 22, 45]:
        label += '_' + hashlib.sha256(json.dumps([args.region, args.input_region]).encode()).hexdigest()[:8]
    args.output_dir = (args.output_dir or ROOT / 'GNSS_analysis/output/run3_gpsgridder_eu' / label).resolve()
    return data, metadata


def main():
    args = arguments()
    data, block = load_inputs(args)
    total = 2 * len(data)
    print(f'Run 3: {len(data)} input representatives; {total} total modes.', flush=True)
    if args.show_mode_counts:
        for percent in (100, 75, 50, 25):
            print(f'-Cn{percent}% -> {requested_mode_count(len(data), percent)} requested modes')
        return
    config = configuration(args, data, block)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    config_path = output / 'run_config.json'
    gmt, env = locate_gmt(args.gmt), gmt_environment()
    command = [gmt, 'gpsgridder', str(output / 'gnss_input_eu.txt'),
               '-R' + '/'.join(f'{v:g}' for v in args.region), f'-I{args.spacing}', '-fg', '-W', '-r',
               f'-Cn{args.retain_percent:g}%']
    if args.poisson is not None:
        command.append(f'-S{args.poisson:g}')
    if args.fudge_km is not None:
        command.append(f'-Fd{args.fudge_km:g}')
    command += [f'-G{output / "velocity_eu.nc"}', f'-E{output / "gmt_station_misfit.txt"}', '-V']
    content = io.StringIO()
    np.savetxt(content, data, fmt='%.8f', header='lon lat VE_mm_yr VN_mm_yr sigmaE_mm_yr sigmaN_mm_yr')
    expected_input = content.getvalue()
    completed = all((output / name).is_file() for name in
                    ('run_config.json', 'gnss_input_eu.txt', 'velocity_eu_u.nc', 'velocity_eu_v.nc', 'gmt_station_misfit.txt', 'gpsgridder.log'))
    matches = completed and json.loads(config_path.read_text())['parameters'] == config and (output / 'gnss_input_eu.txt').read_text() == expected_input
    if args.reuse_grids and not matches:
        raise ValueError('No existing fit matches the current input and settings')
    if completed and not matches and not args.rerun:
        raise ValueError('Existing output differs from this configuration; choose another output directory or use --rerun')
    if matches and not args.rerun:
        print('Reusing identical completed gpsgridder fit.', flush=True)
        if args.verbose:
            print((output / 'gpsgridder.log').read_text(), flush=True)
    else:
        (output / 'gnss_input_eu.txt').write_text(expected_input)
        version = run_command([gmt, '--version'], env, capture_output=True, text=True).stdout.strip()
        config_path.write_text(json.dumps({'stage': 3, 'script': Path(__file__).name,
            'parameters': config, 'input_dir': str(args.input_dir),
            'gmt_version': version, 'command': shlex.join(command)}, indent=2) + '\n')
        print(shlex.join(command), flush=True)
        started = time.monotonic()
        run_gpsgridder(command, env, output / 'gpsgridder.log', verbose=args.verbose)
        print(f'GMT fit completed in {time.monotonic()-started:.1f} s.', flush=True)
    report_used_modes(output, len(data), args.retain_percent)
    print(json.dumps(diagnostics(gmt, output, args, data, env), indent=2), flush=True)
    print(f'Output: {output}', flush=True)


if __name__ == '__main__':
    main()
