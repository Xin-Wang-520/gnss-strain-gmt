#!/usr/bin/env python3
"""Run 5: replot saved run4 fields with bilinear display interpolation only.

python3 run5_plot_strain_eu.py
Writes into display_bilinear/ under the selected run4 directory. No GMT
calculation, numerical filtering, strain recomputation or grid modification.
"""
import argparse
import json
from pathlib import Path
from run4_calculate_strain_eu import ROOT, plot_results, sha256
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / 'GNSS_analysis/output/run4_strain_eu/bm30m_c25_i20m_pixel')
    parser.add_argument('--display-factor', type=int, default=10)
    parser.add_argument('--dpi', type=int, default=220)
    args = parser.parse_args()
    if not 1 <= args.display_factor <= 20 or args.dpi <= 0:
        parser.error('Use display factor 1..20 and positive dpi')
    source = args.input_dir.resolve()
    manifest = json.loads((source / 'run4_summary.json').read_text())
    grid_file = source / 'strain_eu.npz'
    originals = [grid_file, source / 'run4_summary.json', *source.glob('*_eu.nc'),
                 *source.glob('run4_*.png'), *source.glob('run4_*.pdf'), *source.glob('run4_*.svg')]
    before = {str(path): sha256(path) for path in originals}
    with np.load(grid_file) as data:
        lon, lat = data['longitude'], data['latitude']
        fields = {k: data[k] for k in data.files if k not in ('longitude', 'latitude', 'edge_boundary')}
    output = source / 'display_bilinear'
    output.mkdir(exist_ok=True)
    percentile = manifest['plot_colors']['percentile']
    colors = plot_results(output, lon, lat, fields, manifest['run3_parameters']['region'],
                          manifest['earth_radius_km_for_distance_scale'], percentile,
                          args.dpi, display_factor=args.display_factor)
    assert colors == manifest['plot_colors'], 'Native color limits changed'
    assert before == {str(path): sha256(path) for path in originals}, 'An original file changed'
    record = {'source_npz': str(grid_file), 'source_sha256': before[str(grid_file)],
              'native_grid_arcmin': manifest['run3_parameters']['spacing_degree'] * 60,
              'native_grid_shape_lat_lon': [len(lat), len(lon)],
              'display_factor': args.display_factor,
              'interpolation': 'separable bilinear in lon/lat; nearest value within outer half-cells',
              'display_only': True, 'original_files_unchanged_verified': True,
              'color_limits_from_native_fields': colors,
              'note': 'Display interpolation does not add observational information or physical resolution.'}
    (output / 'display_plot_summary.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'Plots saved: {output}\nVerified: native grids, original figures and color limits unchanged.')


if __name__ == '__main__':
    main()
