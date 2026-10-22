#!/usr/bin/env python3
"""Run 4: Xu GMT strain workflow on the accepted Run 3 Eurasia-fixed grids.

Usage: python3 run4_calculate_strain_eu.py
Input: pixel-registered Run 3 VE/VN NetCDF grids, mm/yr.
The calculation commands follow xu/plot_tensor.csh, plot_shear.csh and
plot_dilatation.csh: convert velocity to m/yr, grdgradient -A270/-A180 -fg,
symmetrize the two cross-gradients, then compute principal strain rates.
Maximum shear = abs(e1-e2)/2; dilatation = e1+e2.

No spherical curvature terms, additional velocity smoothing or terrain mask.
The default R=6378 km sets GMT's metres-per-degree distance scale only.
Output units: nstrain/yr; EN is tensor shear, engineering shear is 2 EN.
The second_invariant name follows Xu/RSE and denotes the Frobenius norm.
Python orchestrates GMT and makes figures; GMT calculates the strain fields.
Historical scripts and tests are retained under GNSS_analysis/archive/.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import shutil
import shlex

from run1_plot_gnss import ROOT, draw_context, mercator_y
from GNSS_analysis.gnss_processing import locate_gmt, gmt_environment
import numpy as np
import matplotlib
matplotlib.use("Agg")

DEFAULT_INPUT = ROOT / "GNSS_analysis/output/run3_gpsgridder_eu/bm30m_c25_i20m_pixel"
RATE_KEYS = ("ee", "nn", "en", "max_shear", "dilatation", "second_invariant")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--earth-radius-km", type=float, default=6378)
    parser.add_argument("--plot-percentile", type=float, default=99,
                        help="Shared color limits; saved data are never clipped")
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--gmt", default=None)
    args = parser.parse_args()
    if args.earth_radius_km <= 0 or not np.isfinite(args.earth_radius_km):
        parser.error("Earth radius must be finite and positive")
    if not 90 <= args.plot_percentile <= 100 or args.dpi <= 0:
        parser.error("Use 90 <= plot-percentile <= 100 and positive dpi")
    args.input_dir = args.input_dir.resolve()
    if args.output_dir is None:
        name = args.input_dir.name
        if args.earth_radius_km != 6378:
            name += "_R" + f"{args.earth_radius_km:g}".replace(".", "p")
        args.output_dir = ROOT / "GNSS_analysis/output/run4_strain_eu" / name
    args.output_dir = args.output_dir.resolve()
    return args


def read_grid_binary(gmt, path, env):
    """Export native GMT grid values without text rounding."""
    data = np.frombuffer(subprocess.check_output(
        [gmt, "grd2xyz", str(path), "-bo3d"], env=env), dtype=np.float64).reshape(-1, 3)
    longitude = np.unique(data[:, 0])
    latitude = data[::len(longitude), 1]
    return longitude, latitude, data[:, 2].reshape(len(latitude), len(longitude))


def xu_strain(gmt, east_path, north_path, workdir, radius_km, env):
    """Execute Xu's GMT arithmetic, keeping source-compatible intermediate names.

    exy.grd and eyx.grd are *velocity gradients*. e12.grd is the tensor EN.
    The original terrain-mask and plotting sections are omitted: the accepted
    GNSS result is currently unmasked. This function does not run gpsgridder.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    commands = []
    logs = []

    def call(*arguments):
        cmd = [gmt, *map(str, arguments)]
        commands.append(shlex.join(cmd))
        result = subprocess.run(cmd, cwd=workdir, env=env, text=True, capture_output=True)
        logs.append("$ " + shlex.join(cmd) + "\n" + result.stdout + result.stderr)
        (workdir / "gmt_execution.log").write_text("\n".join(logs))
        if result.returncode:
            raise RuntimeError(f"GMT failed: {shlex.join(cmd)}\n{result.stderr}")

    call("set", "PROJ_ELLIPSOID", f"{radius_km*1000:.9g}", "PROJ_MEAN_RADIUS", "authalic")
    settings = subprocess.check_output([gmt, "get", "PROJ_ELLIPSOID", "PROJ_MEAN_RADIUS"],
                                       cwd=workdir, env=env, text=True).strip()
    call("grdmath", east_path, "1000", "DIV", "=", "GNSS_um.grd")
    call("grdmath", north_path, "1000", "DIV", "=", "GNSS_vm.grd")
    for name in ("GNSS_um.grd", "GNSS_vm.grd"):
        call("grdedit", name, "-D+dm/yr+tEurasia-fixed GNSS velocity")
    # Same four directional derivatives as the provided Xu scripts.
    call("grdgradient", "GNSS_um.grd", "-A270", "-fg", "-Gexx.grd")
    call("grdgradient", "GNSS_um.grd", "-A180", "-fg", "-Gexy.grd")
    call("grdgradient", "GNSS_vm.grd", "-A270", "-fg", "-Geyx.grd")
    call("grdgradient", "GNSS_vm.grd", "-A180", "-fg", "-Geyy.grd")
    # plot_tensor.csh evaluates this norm before converting components to nstrain.
    call("grdmath", "exy.grd", "eyx.grd", "ADD", "2", "DIV", "=", "e12.grd")
    call("grdmath", "exx.grd", "exx.grd", "MUL", "eyy.grd", "eyy.grd", "MUL", "ADD",
         "e12.grd", "e12.grd", "MUL", "2", "MUL", "ADD", "SQRT", "1E9", "MUL", "=", "e_2nd.grd")
    for name in ("exx.grd", "exy.grd", "eyx.grd", "eyy.grd"):
        call("grdmath", name, "1E9", "MUL", "=", name)
    call("grdmath", "exy.grd", "eyx.grd", "ADD", "2", "DIV", "=", "e12.grd")
    call("grdmath", "2", "e12.grd", "MUL", "exx.grd", "eyy.grd", "SUB", "ATAN2", "2", "DIV", "=", "a.grd")
    common = ["exx.grd", "eyy.grd", "SUB", "2", "POW", "0.25", "MUL", "e12.grd", "2", "POW", "ADD", "SQRT"]
    call("grdmath", *common, "exx.grd", "eyy.grd", "ADD", "2", "DIV", "ADD", "=", "e1.grd")
    call("grdmath", *common, "NEG", "exx.grd", "eyy.grd", "ADD", "2", "DIV", "ADD", "=", "e2.grd")
    call("grdmath", "e1.grd", "e2.grd", "SUB", "ABS", "2", "DIV", "=", "shear.grd")
    call("grdmath", "e1.grd", "e2.grd", "ADD", "=", "dilatation.grd")
    (workdir / "run4_xu_commands.sh").write_text(
        "#!/usr/bin/env bash\nset -euo pipefail\ncd -- " + shlex.quote(str(workdir)) + "\n" + "\n".join(commands) + "\n")
    grid_names = {"ee": "exx.grd", "nn": "eyy.grd", "en": "e12.grd", "max_shear": "shear.grd",
                  "dilatation": "dilatation.grd", "second_invariant": "e_2nd.grd",
                  "principal_max": "e1.grd", "principal_min": "e2.grd",
                  "principal_angle_from_east_rad": "a.grd"}
    fields = {}
    for key, name in grid_names.items():
        lon, lat, fields[key] = read_grid_binary(gmt, workdir / name, env)
    a, b = fields["en"]**2, ((fields["ee"]-fields["nn"])/2)**2
    fields.update(shear_en_squared=a, shear_normal_difference_squared=b,
                  normal_difference_fraction=np.divide(b, a+b, out=np.full_like(b, np.nan), where=a+b>0))
    return lon, lat, fields, grid_names, settings


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def field_statistics(values):
    finite = values[np.isfinite(values)]
    if not len(finite):
        return {"finite_count": 0}
    return {"finite_count": int(len(finite)), "min": float(finite.min()),
            "p01": float(np.percentile(finite, 1)), "median": float(np.median(finite)),
            "p99": float(np.percentile(finite, 99)), "max": float(finite.max())}


def bilinear_display_grid(longitude, latitude, values, factor):
    """Interpolate pixel values for display only; constant at outer half-cells.

    Interpolation is in geographic coordinates before Mercator plotting.
    The input arrays are never modified. No new physical resolution is implied.
    """
    if factor == 1:
        return longitude, latitude, values
    if factor < 1 or int(factor) != factor:
        raise ValueError("Display factor must be a positive integer")
    dx = float(np.median(np.diff(longitude)))
    dy = float(np.median(abs(np.diff(latitude))))
    if dx <= 0 or dy <= 0 or not np.allclose(np.diff(longitude), dx) or not np.allclose(abs(np.diff(latitude)), dy):
        raise ValueError("Display interpolation expects a regular geographic grid")
    nx, ny = len(longitude) * factor, len(latitude) * factor
    fine_lon = longitude[0] - dx / 2 + (np.arange(nx) + .5) * dx / factor
    fine_lat = latitude.min() - dy / 2 + (np.arange(ny) + .5) * dy / factor
    if latitude[0] > latitude[-1]:
        coarse_lat, coarse_values = latitude[::-1], values[::-1]
        fine_lat = fine_lat[::-1]
    else:
        coarse_lat, coarse_values = latitude, values
    along_lon = np.stack([np.interp(fine_lon, longitude, row) for row in coarse_values])
    fine_values = np.stack([np.interp(fine_lat, coarse_lat, column) for column in along_lon.T], axis=1)
    return fine_lon, fine_lat, fine_values


def plot_results(output, longitude, latitude, fields, region, radius_km, percentile, dpi,
                 display_factor=1):
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    display_note = "No added smoothing"
    if display_factor > 1:
        native_arcmin = abs(float(longitude[1] - longitude[0])) * 60
        display_note = f"Native grid {native_arcmin:g} arcmin | Bilinear display interpolation only"

    def draw_field(ax, values, norm, cmap):
        plot_lon, plot_lat, plot_values = bilinear_display_grid(longitude, latitude, values, display_factor)
        lon, lat = np.meshgrid(plot_lon, plot_lat)
        return ax.pcolormesh(lon, mercator_y(lat), plot_values, shading="auto", rasterized=True,
                             norm=norm, cmap=cmap, zorder=0.5)
    signed = np.concatenate([fields[key].ravel() for key in ("ee", "nn", "en", "dilatation")])
    positive = np.concatenate([fields[key].ravel() for key in ("max_shear", "second_invariant")])
    signed_limit = max(1.0, float(np.nanpercentile(np.abs(signed), percentile)))
    positive_limit = max(1.0, float(np.nanpercentile(positive, percentile)))
    signed_norm = TwoSlopeNorm(vmin=-signed_limit, vcenter=0, vmax=signed_limit)
    positive_norm = Normalize(vmin=0, vmax=positive_limit)
    fig, axes = plt.subplots(2, 3, figsize=(20, 10), layout="constrained")
    fig.suptitle(f"Run 4 | GNSS horizontal strain rates | Eurasia-fixed\n"
                 f"Xu GMT workflow | R={radius_km:g} km for distance scale | {display_note}", fontsize=16)
    specs = [("ee", "a  EE: east-west extension / shortening", False),
             ("nn", "b  NN: north-south extension / shortening", False),
             ("en", "c  EN: tensor shear (half engineering shear)", False),
             ("max_shear", "d  Maximum shear: half principal-rate difference", True),
             ("dilatation", "e  Horizontal dilatation: EE + NN", False),
             ("second_invariant", "f  Tensor magnitude: sqrt(EE^2 + NN^2 + 2 EN^2)", True)]
    for ax, (key, title, is_positive) in zip(axes.ravel(), specs):
        draw_context(ax, region, ROOT / "GNSS_analysis/data/basemap")
        im = draw_field(ax, fields[key], positive_norm if is_positive else signed_norm,
                        "YlOrRd" if is_positive else "RdBu_r")
        ax.set_title(title, loc="left", fontsize=10)
        fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8, extend="max" if is_positive else "both", label="nstrain/yr")
    footer = f"Positive EE/NN/dilatation = extension. Shared color limits use P{percentile:g}; saved grids are not clipped. Outer-cell gradients depend on GMT boundary conditions."
    if display_factor > 1:
        footer = f"Bilinear interpolation for display only; original grids and physical resolution unchanged. Color limits use native-grid P{percentile:g}."
    fig.supxlabel(footer, fontsize=10)
    fig.savefig(output / "run4_strain_eu.png", dpi=dpi, facecolor="white")
    fig.savefig(output / "run4_strain_eu.pdf", dpi=dpi, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5), layout="constrained")
    contribution_title = "Maximum-shear components | max_shear^2 = A + B"
    if display_factor > 1:
        contribution_title += "\n" + display_note
    fig.suptitle(contribution_title, fontsize=16)
    pieces = [(np.sqrt(fields["shear_en_squared"]), "sqrt(A) = abs(EN)", positive_norm, "YlOrRd", "nstrain/yr"),
              (np.sqrt(fields["shear_normal_difference_squared"]), "sqrt(B) = abs(EE - NN) / 2", positive_norm, "YlOrRd", "nstrain/yr"),
              (fields["normal_difference_fraction"], "Normal-difference fraction: B / (A + B)", Normalize(0, 1), "viridis", "Fraction")]
    for ax, (values, title, norm, cmap, unit) in zip(axes, pieces):
        draw_context(ax, region, ROOT / "GNSS_analysis/data/basemap")
        im = draw_field(ax, values, norm, cmap)
        ax.set_title(title, loc="left", fontsize=11)
        fig.colorbar(im, ax=ax, pad=0.02, shrink=0.8, extend="neither" if unit == "Fraction" else "max", label=unit)
    fig.supxlabel("A large B fraction means the EE-NN difference dominates; inspect EE and NN separately before attributing it to NN. Fractions at very small shear are not diagnostic.", fontsize=10)
    fig.savefig(output / "run4_shear_contributions_eu.png", dpi=dpi, facecolor="white")
    fig.savefig(output / "run4_shear_contributions_eu.pdf", dpi=dpi, facecolor="white")
    plt.close(fig)
    return {"signed_color_limit_nstrain_yr": signed_limit,
            "positive_color_limit_nstrain_yr": positive_limit, "percentile": percentile}


def main():
    args = arguments()
    source = args.input_dir
    config_path = source / "run_config.json"
    config = json.loads(config_path.read_text())
    parameters = config["parameters"]
    if parameters["frame"] != "ITRF2014 Eurasia-fixed":
        raise ValueError("Run 4 EU requires Eurasia-fixed velocities")
    if parameters.get("registration") != "pixel (-r)":
        raise ValueError("Expected accepted pixel-registered Run 3 grids")
    region, increment = parameters["region"], parameters["gmt_increment"]
    gmt, env = locate_gmt(args.gmt), gmt_environment()
    paths = [source / "velocity_eu_u.nc", source / "velocity_eu_v.nc"]
    longitude, latitude, ve = read_grid_binary(gmt, paths[0], env)
    lon2, lat2, vn = read_grid_binary(gmt, paths[1], env)
    if not (np.array_equal(longitude, lon2) and np.array_equal(latitude, lat2)):
        raise ValueError("East/north velocity grid axes differ")
    # Verify pixel centers and spacing against the recorded interpolation region.
    step = parameters["spacing_degree"]
    if not (np.allclose(np.diff(longitude), step) and np.allclose(np.abs(np.diff(latitude)), step)):
        raise ValueError("Unexpected grid spacing")
    if not (np.allclose([longitude.min()-step/2, longitude.max()+step/2], region[:2]) and
            np.allclose([latitude.min()-step/2, latitude.max()+step/2], region[2:])):
        raise ValueError("Grid is not pixel-centered on the recorded region")
    print(f"Input: {source}\nGrid: {len(longitude)} x {len(latitude)}; running Xu GMT strain workflow.", flush=True)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    work = output / "xu_intermediate"
    xlon, xlat, fields, grids, settings = xu_strain(
        gmt, paths[0], paths[1], work, args.earth_radius_km, env)
    if not (np.array_equal(longitude, xlon) and np.array_equal(latitude, xlat)):
        raise ValueError("GMT output axes differ from input pixel centers")
    for key in RATE_KEYS:
        if not np.isfinite(fields[key]).all():
            raise ValueError(f"Non-finite strain field: {key}")
    edge = np.zeros(ve.shape, dtype=bool)
    edge[[0, -1], :] = True
    edge[:, [0, -1]] = True
    np.savez_compressed(output / "strain_eu.npz", longitude=longitude, latitude=latitude,
                        edge_boundary=edge, **fields)
    for key, name in grids.items():
        destination = output / f"{key}_eu.nc"
        shutil.copy2(work / name, destination)
        units = "radian" if key == "principal_angle_from_east_rad" else "nstrain/yr"
        subprocess.run([gmt, "grdedit", str(destination), "-Cp",
            f"-D+d{units}+t{key.replace('_',' ')} (Xu GMT method)+rXu directional gradients and principal-rate arithmetic; planar metric"],
            env=env, check=True, capture_output=True, text=True)
        print(f"Saved {destination.name}", flush=True)
    # This fraction is an additional diagnostic for the boss's NN/EE question.
    subprocess.run([gmt, "grdmath", str(output / "ee_eu.nc"), str(output / "nn_eu.nc"),
        "SUB", "2", "POW", "4", "DIV", "=", str(work / "normal_difference_squared.grd")],
        env=env, check=True, capture_output=True, text=True)
    subprocess.run([gmt, "grdmath", str(output / "en_eu.nc"), "2", "POW",
        str(work / "normal_difference_squared.grd"), "ADD", "0", "NAN", "=", str(work / "shear_squared.grd")],
        env=env, check=True, capture_output=True, text=True)
    subprocess.run([gmt, "grdmath", str(work / "normal_difference_squared.grd"),
        str(work / "shear_squared.grd"), "DIV", "=", str(output / "normal_difference_fraction_eu.nc")],
        env=env, check=True, capture_output=True, text=True)
    subprocess.run([gmt, "grdedit", str(output / "normal_difference_fraction_eu.nc"),
        "-D+d1+tNormal difference fraction+rNormal difference squared contribution divided by total shear squared"],
        env=env, check=True, capture_output=True, text=True)
    colors = plot_results(output, longitude, latitude, fields, region,
                           args.earth_radius_km, args.plot_percentile, args.dpi)
    statistics = {key: field_statistics(fields[key]) for key in RATE_KEYS + ("normal_difference_fraction",)}
    source_scripts = [ROOT / "xu" / name for name in ("plot_tensor.csh", "plot_shear.csh", "plot_dilatation.csh")]
    manifest = {"script": Path(__file__).name, "method": "Xu GMT directional gradients and principal strain rates",
                "input_dir": str(source),
                "input_files": [{"path": str(path), "sha256": sha256(path)} for path in paths],
                "source_scripts": [{"path": str(path), "sha256": sha256(path)} for path in source_scripts],
                "run3_parameters": parameters, "frame": parameters["frame"],
                "coordinates": "East-North; geographic grids with GMT flat-distance metric",
                "earth_radius_km_for_distance_scale": args.earth_radius_km,
                "gmt_settings": settings, "gmt_version": subprocess.check_output([gmt, "--version"], text=True).strip(),
                "input_velocity_units": "mm/yr", "intermediate_velocity_units": "m/yr",
                "output_strain_units": "nstrain/yr = 1e-9 / yr",
                "shear_convention": "EN=e12=(du/dy+dv/dx)/2; exy.grd is du/dy, not EN",
                "max_shear_definition": "abs(e1-e2)/2, calculated by GMT as in Xu plot_shear.csh",
                "dilatation_definition": "e1+e2, calculated by GMT as in Xu plot_dilatation.csh",
                "second_invariant_definition": "sqrt(EE^2+NN^2+2EN^2), Frobenius norm as in Xu plot_tensor.csh",
                "spherical_curvature_terms": False, "added_smoothing": False,
                "observation_distance_mask": False, "terrain_mask": False,
                "derivatives": "GMT grdgradient -A270/-A180 -fg; default natural boundary conditions",
                "grid_shape_lat_lon": list(ve.shape), "grid_registration": parameters["registration"],
                "principal_angle": "a.grd = 0.5 atan2(2EN, EE-NN), radians counterclockwise from east",
                "command_log": str(work / "run4_xu_commands.sh"), "plot_colors": colors,
                "adaptations": ["Use current Run 3 EU VE/VN, converted from mm/yr to m/yr",
                                "Set metric radius explicitly to 6378 km by default, matching previous RSE scale",
                                "No terrain mask: retain the currently accepted complete GNSS region",
                                "Replace California plotting overlays with current study-region maps",
                                "Add shear-contribution diagnostics; core rates are computed by GMT"],
                "limitations": ["No new velocity smoothing: local noise can amplify in derivatives",
                                "Sparse-station and boundary areas remain unmasked",
                                "No propagation of GNSS/interpolation uncertainty or independent validation"],
                "fields": statistics}
    (output / "run4_summary.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({key: statistics[key] for key in ("nn", "max_shear", "dilatation")}, indent=2), flush=True)
    print(f"Results: {output}", flush=True)


if __name__ == "__main__":
    main()
