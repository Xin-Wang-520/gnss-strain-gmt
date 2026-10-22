"""Parameter sweeps for the three explicit GNSS processing stages."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time
from run1_plot_gnss import ROOT, draw_context, mercator_y
from .gnss_processing import token
from run4_calculate_strain_eu import RATE_KEYS
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm


def field_statistics(values):
    return {'median': float(np.median(values)), 'p99_abs': float(np.percentile(abs(values), 99)),
            'min': float(values.min()), 'max': float(values.max())}


def compare(cases, output, kind, dpi):
    records = {}
    common_reference = None
    for case in cases:
        folder = Path(case['folder'])
        config = json.loads((folder / 'run3_gpsgridder/run_config.json').read_text())
        pars = config['parameters']
        ignored = {'retain_percent'} if kind == 'percentages' else {'gmt_increment', 'spacing_degree'}
        if kind == 'blockmedian':
            ignored = {'blockmedian', 'input_count', 'input_values_sha256'}
        common = {k: v for k, v in pars.items() if k not in ignored}
        common['raw_source_sha256'] = pars['blockmedian']['raw_values_sha256']
        if common_reference is None:
            common_reference = common
        if common != common_reference:
            raise ValueError('Parameters other than the intended test variable changed')
        summary = json.loads((folder / 'run4_strain/run4_summary.json').read_text())
        if summary['run3_parameters'] != pars or 'Xu' not in summary['method']:
            raise ValueError('Strain output does not match this velocity configuration')
        with np.load(folder / 'run4_strain/strain_eu.npz') as data:
            lon, lat = data['longitude'], data['latitude']
            fields = {k: data[k] for k in RATE_KEYS}
        x, y = np.meshgrid(lon, lat)
        w, e, s, n = pars['region']
        mask = (x >= w+.5) & (x <= e-.5) & (y >= s+.5) & (y <= n-.5)
        assert mask.any() and all(np.isfinite(v).all() for v in fields.values())
        case.update(lon=lon, lat=lat, fields=fields, mask=mask)
        records[case['name']] = {
            'parameters': pars, 'grid_shape_lat_lon': [len(lat), len(lon)],
            'strain': {k: field_statistics(v[mask]) for k, v in fields.items()},
            'mode_count': json.loads((folder / 'run3_gpsgridder/mode_summary.json').read_text()),
            'fit': json.loads((folder / 'run3_gpsgridder/fit_summary.json').read_text())}
    labels = {'ee': 'EE', 'nn': 'NN', 'en': 'EN: tensor shear', 'max_shear': 'Maximum shear',
              'dilatation': 'Dilatation', 'second_invariant': 'Tensor magnitude'}
    for key in RATE_KEYS:
        positive = key in ('max_shear', 'second_invariant')
        limit = max(1., float(np.percentile(abs(cases[0]['fields'][key][cases[0]['mask']]), 99)))
        norm = Normalize(0, limit) if positive else TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        cols = min(3, len(cases)); rows = (len(cases)+cols-1)//cols
        fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 3.8*rows), squeeze=False, layout='constrained')
        fig.suptitle(f'{labels[key]} | Test {kind} | Xu GMT strain, fixed color scale', fontsize=14)
        for ax, case in zip(axes.ravel(), cases):
            draw_context(ax, (w, e, s, n), ROOT / 'GNSS_analysis/data/basemap')
            x, y = np.meshgrid(case['lon'], case['lat'])
            im = ax.pcolormesh(x, mercator_y(y), case['fields'][key], shading='auto', norm=norm,
                              cmap='YlOrRd' if positive else 'RdBu_r', rasterized=True, zorder=.5)
            ax.set_title(case['name'], fontsize=10)
        for ax in axes.ravel()[len(cases):]:
            ax.set_visible(False)
        fig.colorbar(im, ax=axes.ravel().tolist(), label='nstrain/yr', shrink=.8, pad=.01,
                     extend='max' if positive else 'both')
        fig.supxlabel('Native-grid statistics exclude 0.5 degrees at map edges. In-sample fit and map smoothness are not independent validation.', fontsize=9)
        fig.savefig(output / f'{key}_comparison.png', dpi=dpi, facecolor='white')
        plt.close(fig)
    (output / 'comparison.json').write_text(json.dumps({'test_variable': kind, 'cases': records,
        'notes': ['Only the stated parameter is varied.',
                  'Spacing tests sample different native grid locations, so statistics are not pointwise differences.',
                  'Changing blockmedian also changes representative weights and the default Ff physical distance.',
                  'No independent validation or uncertainty propagation is performed.']}, indent=2)+'\n')
    print(json.dumps({k: v['strain']['max_shear'] for k, v in records.items()}, indent=2), flush=True)


def main(kind):
    defaults = {'blockmedian': ['5', '10', '15', '20', '30'],
                'spacing': ['5', '10', '15', '20'], 'percentages': ['100', '75', '50', '25']}
    parser = argparse.ArgumentParser(description=f'Test {kind}; run2 blockmedian -> run3 gpsgridder -> run4 Xu strain.')
    parser.add_argument('--values', nargs='+', default=defaults[kind], help='Arcminutes for widths/spacings, retained percentages for percentage tests')
    parser.add_argument('--block-spacing', default='30m')
    parser.add_argument('--spacing', default='20m')
    parser.add_argument('--retain-percent', type=float, default=25)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--rerun', action='store_true')
    parser.add_argument('--compare-only', action='store_true')
    parser.add_argument('--dpi', type=int, default=150)
    args = parser.parse_args()
    try:
        numbers = [float(v.rstrip('m')) for v in args.values]
    except ValueError:
        parser.error('Use numeric values, e.g. --values 5 10 15 20')
    if any(not np.isfinite(v) or v <= 0 or (kind == 'percentages' and v > 100) for v in numbers) or args.workers not in (1, 2):
        parser.error('Use positive values (percentages <=100) and one or two workers')
    if len(set(numbers)) != len(numbers):
        parser.error('Duplicate test values are not allowed')
    signature = f'bm{args.block_spacing}_c{token(args.retain_percent)}_i{args.spacing}'
    output = (args.output_dir or ROOT / 'GNSS_analysis/output/tests' / f'{kind}_{signature}').resolve()
    output.mkdir(parents=True, exist_ok=True)
    lock = (output / 'batch.lock').open('w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('Another test batch is already using this output directory')
    cases = []
    for value in numbers:
        block = f'{value:g}m' if kind == 'blockmedian' else args.block_spacing
        spacing = f'{value:g}m' if kind == 'spacing' else args.spacing
        percent = value if kind == 'percentages' else args.retain_percent
        name = f'bm{block}_c{token(percent)}_i{spacing}'
        cases.append({'name': name, 'folder': str(output / name), 'block': block, 'spacing': spacing, 'percent': percent})
    started = time.monotonic()
    stopping = threading.Event(); active = {}; active_lock = threading.Lock()
    def execute(command, log):
        with active_lock:
            if stopping.is_set():
                return
            process = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
            active[process.pid] = process
        try:
            code = process.wait()
            if code and not stopping.is_set():
                raise subprocess.CalledProcessError(code, command)
        finally:
            with active_lock:
                active.pop(process.pid, None)
    def job(case):
        folder = Path(case['folder']); folder.mkdir(parents=True, exist_ok=True)
        block, gps, strain = [folder / name for name in ('run2_blockmedian', 'run3_gpsgridder', 'run4_strain')]
        commands = [
            [sys.executable, str(ROOT / 'run2_blockmedian_eu.py'), '--spacing', case['block'], '--output-dir', str(block)],
            [sys.executable, str(ROOT / 'run3_gpsgridder_eu.py'), '--input-dir', str(block), '--output-dir', str(gps),
             '--spacing', case['spacing'], '--retain-percent', str(case['percent']), '-v'],
            [sys.executable, str(ROOT / 'run4_calculate_strain_eu.py'), '--input-dir', str(gps), '--output-dir', str(strain)]]
        if args.rerun:
            commands[1].append('--rerun')
        (folder / 'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
        print('START ' + case['name'], flush=True)
        with (folder / 'test.log').open('w') as log:
            for command in commands:
                execute(command, log)
        if not stopping.is_set():
            print(f'DONE {case["name"]} ({time.monotonic()-started:.0f}s)', flush=True)
    errors = []
    if not args.compare_only:
        pool = ThreadPoolExecutor(max_workers=args.workers)
        try:
            futures = {pool.submit(job, case): case for case in cases}
            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as exc:
                    errors.append({'case': futures[future]['name'], 'error': str(exc)})
        except KeyboardInterrupt:
            stopping.set()
            with active_lock:
                for process in active.values():
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
            pool.shutdown(wait=True, cancel_futures=True)
            raise SystemExit('Stopped; completed results retained')
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
    (output / 'batch_summary.json').write_text(json.dumps({'cases': cases, 'errors': errors,
        'elapsed_seconds': time.monotonic()-started}, indent=2)+'\n')
    if errors:
        raise SystemExit(f'Test cases failed; inspect logs under {output}')
    compare(cases, output, kind, args.dpi)
    print(f'Test results: {output}', flush=True)

