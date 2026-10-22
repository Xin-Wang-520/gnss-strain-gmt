"""Shared GMT preprocessing, velocity diagnostics and plotting helpers.

User entry points are run2_blockmedian_eu.py, run3_gpsgridder_eu.py and
run4_calculate_strain_eu.py. Numerical methods are preserved from the
validated original implementation.
"""
from pathlib import Path
import hashlib
import io
import json
import os
import re
import shlex
import shutil
import subprocess
from run1_plot_gnss import ROOT, draw_context, mercator_y
import numpy as np
import matplotlib
matplotlib.use("Agg")

def increment_degrees(text):
    """Convert the angular increments supported by this script to degrees."""
    text = str(text).strip()
    factor = {"d": 1.0, "m": 1/60, "s": 1/3600}
    if text[-1:] in factor:
        return float(text[:-1]) * factor[text[-1]]
    return float(text)


def token(value):
    return f"{value:g}".replace("-", "m").replace(".", "p")


def configuration(args, data, block_metadata=None):
    result = {"frame": "ITRF2014 Eurasia-fixed", "version": "AHB v6.4", "region": args.region,
            "method": "Run 6.1 gpsgridder settings", "registration": "pixel (-r)",
            "input_region": args.input_region, "spacing_degree": args.spacing_degree,
            "gmt_increment": args.spacing,
            "poisson_2d": args.poisson, "retain_percent": args.retain_percent,
            "poisson_note": "None means no -S option: use GMT default, 0.5 in GMT 6.4",
            "fudge_km": args.fudge_km, "weights": "-W: standard-deviation input, GMT default sigma weighting",
            "fudge_note": "None means no -F option: use GMT default, -Ff0.01 in GMT 6.4",
            "input_count": len(data),
            "input_values_sha256": hashlib.sha256(np.asarray(data, dtype="<f8").tobytes()).hexdigest(),
            "distance_mode": "GMT -fg: single-region flat-Earth approximation"}
    if block_metadata is not None:
        result["blockmedian"] = block_metadata
        result["input_kind"] = "block representative velocities, not individual stations"
    return result


def blockmedian_input(gmt, data, region, increment, output, env):
    """Component-wise unweighted medians in identical geographic blocks.

    GMT's default median x/y depends only on the shared input locations, so
    VE, VN and sigma runs retain identical representative coordinates. Do not
    use -Q, which would select different positions for different components.
    Sigma medians retain a typical uncertainty scale; they are not a formal
    propagated uncertainty of the median or reduced by sqrt(block count).
    """
    output.mkdir(parents=True, exist_ok=True)
    original = output / "gnss_before_blockmedian_eu.txt"
    np.savetxt(original, data, fmt="%.12g", header="lon lat VE VN sigmaE sigmaN; velocities and sigma in mm/yr")
    tables, commands = [], []
    counts = None
    common = ["-R"+"/".join(f"{v:g}" for v in region), f"-I{increment}", "-fg", "-r"]
    for column, name in ((2, "ve"), (3, "vn"), (4, "sigma_e"), (5, "sigma_n")):
        cmd = [gmt, "blockmedian", str(original), *common, f"-i0,1,{column}"]
        if column == 2:
            cmd += ["-Wo", "-bo4d"]
            ncol = 4
        else:
            cmd += ["-bo3d"]
            ncol = 3
        result = run_command(cmd, env, capture_output=True)
        table = np.frombuffer(result.stdout, dtype=np.float64).reshape(-1, ncol).copy()
        if len(table) < 1:
            raise ValueError("blockmedian produced no constraints")
        if tables:
            np.testing.assert_array_equal(table[:, :2], tables[0][:, :2])
        else:
            counts = table[:, 3].astype(int)
        tables.append(table)
        commands.append(shlex.join(cmd))
        np.savetxt(output / f"blockmedian_{name}.txt", table, fmt="%.12g",
                   header="lon lat median_value"+(" count" if column == 2 else ""))
    reduced = np.column_stack((tables[0][:, :2], *[t[:, 2] for t in tables]))
    if (not np.isfinite(reduced).all() or np.any(reduced[:, 4:6] <= 0) or
            counts.sum() != len(data) or len(np.unique(reduced[:, :2], axis=0)) != len(reduced)):
        raise ValueError("Invalid block medians, nonpositive sigma, lost observations or duplicate block positions")
    metadata = {"gmt_increment": increment, "region": region, "registration": "pixel (-r)",
        "raw_station_count": len(data), "occupied_block_count": len(reduced),
        "merged_block_count": int((counts > 1).sum()), "largest_block_count": int(counts.max()),
        "positions": "unweighted median longitude and latitude in each block, GMT default; not -Q",
        "velocities": "separate unweighted medians of VE and VN over exactly the same blocks",
        "sigma": "component-wise median of station sigma; representative scale, not propagated median uncertainty",
        "raw_values_sha256": hashlib.sha256(np.asarray(data, dtype="<f8").tobytes()).hexdigest()}
    np.savetxt(output / "blockmedian_counts.txt", np.column_stack((reduced[:, :2], counts)),
               fmt="%.12g", header="representative_lon representative_lat original_station_count")
    (output / "blockmedian_summary.json").write_text(json.dumps(metadata, indent=2)+"\n")
    (output / "blockmedian_commands.txt").write_text("\n".join(commands)+"\n")
    return reduced, metadata


def locate_gmt(requested):
    executable = requested or shutil.which("gmt")
    fallback = Path("/Applications/GMT-6.4.0.app/Contents/Resources/bin/gmt")
    if not executable and fallback.is_file():
        executable = str(fallback)
    if not executable:
        raise RuntimeError("GMT 6.x is required; use --gmt /absolute/path/to/gmt")
    return str(Path(executable).resolve())


def gmt_environment():
    env = os.environ.copy()
    cache = ROOT / "GNSS_analysis" / ".cache" / "gmt_processing"
    cache.mkdir(parents=True, exist_ok=True)
    env["GMT_USERDIR"] = str(cache)
    env.setdefault("VECLIB_MAXIMUM_THREADS", "4")
    env.setdefault("OMP_NUM_THREADS", "4")
    return env


def run_command(cmd, env, **kwargs):
    return subprocess.run(cmd, env=env, check=True, **kwargs)


def requested_mode_count(station_count, percent):
    """GMT 6.4 -Cn count: nearest integer, exact halves rounded to even.

    This is a requested truncation count, not a claim of numerical matrix rank.
    After solving, the GMT log supplies the actual count used.
    """
    return int(np.rint((2 * station_count) * (percent / 100)))


def run_gpsgridder(command, env, log_path, verbose=False):
    """Keep the existing log; optionally tee merged GMT output to the terminal."""
    with log_path.open("w") as log:
        if not verbose:
            run_command(command, env, stdout=log, stderr=log)
            return
        process = subprocess.Popen(command, env=env, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, bufsize=1)
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end="", flush=True)
            code = process.wait()
            if code:
                raise subprocess.CalledProcessError(code, command)
        except BaseException:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise
        finally:
            process.stdout.close()


def report_used_modes(output, station_count, percent):
    log_path = output / "gpsgridder.log"
    log_text = log_path.read_text() if log_path.is_file() else ""
    matches = re.findall(r"\[(\d+) of (\d+) eigen-values used\]", log_text)
    record = {"input_point_count": station_count, "total_modes": 2 * station_count,
              "requested_percent": percent,
              "requested_mode_count": requested_mode_count(station_count, percent),
              "actual_modes_used": None, "source": str(log_path)}
    if matches:
        used, total = map(int, matches[-1])
        if total != 2 * station_count:
            raise ValueError("GMT log mode count does not match current input")
        record["actual_modes_used"] = used
        print(f"GMT actual modes used: {used} / {total} ({percent:g}% requested).", flush=True)
    else:
        print("GMT actual mode count was not found in the log; requested count is not treated as verified.", flush=True)
    (output / "mode_summary.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def load_grid(gmt, path, env):
    result = run_command([gmt, "grd2xyz", str(path), "--FORMAT_FLOAT_OUT=%.12g"],
                         env, capture_output=True, text=True)
    xyz = np.loadtxt(io.StringIO(result.stdout))
    longitude = np.unique(xyz[:, 0])
    nlon = len(longitude)
    latitude = xyz[::nlon, 1]
    assert len(xyz) == len(longitude) * len(latitude)
    assert np.allclose(xyz[:nlon, 0], longitude)
    return longitude, latitude, xyz[:, 2].reshape(len(latitude), nlon)


def nearest_distance_km(lon, lat, data):
    """Nearest-station great-circle distance, chunked to limit memory."""
    def unit_xyz(longitude, latitude):
        lam, phi = np.deg2rad(longitude), np.deg2rad(latitude)
        return np.column_stack((np.cos(phi)*np.cos(lam), np.cos(phi)*np.sin(lam), np.sin(phi)))
    stations = unit_xyz(data[:, 0], data[:, 1])
    grid = unit_xyz(lon.ravel(), lat.ravel())
    distance = np.empty(len(grid))
    for start in range(0, len(grid), 512):
        cosine = (grid[start:start+512] @ stations.T).max(axis=1)
        distance[start:start+512] = 6371.0088 * np.arccos(np.clip(cosine, -1, 1))
    return distance.reshape(lon.shape)


def label_grid(gmt, path, title, units, env):
    run_command([gmt, "grdedit", str(path), "-Cp",
                 f"-D+d{units}+t{title}+rAHB v6.4 GNSS; ITRF2014 Eurasia-fixed; velocity interpolation, not strain"],
                 env, capture_output=True, text=True)


def save_grid(gmt, path, lon, lat, values, args, env, units="mm/yr"):
    content = io.StringIO()
    np.savetxt(content, np.column_stack((lon.ravel(), lat.ravel(), values.ravel())), fmt="%.10g")
    run_command([gmt, "xyz2grd", "-R"+"/".join(f"{v:g}" for v in args.region),
                 f"-I{args.spacing}", "-fg", "-r", f"-G{path}"], env, input=content.getvalue(), text=True,
                 capture_output=True)
    label_grid(gmt, path, path.stem, units, env)


def diagnostics(gmt, output, args, data, env):
    # GMT 6.4's default remark calls these "strain components" even when
    # velocities are supplied. Correct the metadata without changing values.
    label_grid(gmt, output / "velocity_eu_u.nc", "Eastward velocity (Eurasia-fixed)", "mm/yr", env)
    label_grid(gmt, output / "velocity_eu_v.nc", "Northward velocity (Eurasia-fixed)", "mm/yr", env)
    longitude, latitude, ve = load_grid(gmt, output / "velocity_eu_u.nc", env)
    lon2, lat2, vn = load_grid(gmt, output / "velocity_eu_v.nc", env)
    assert np.array_equal(longitude, lon2) and np.array_equal(latitude, lat2)
    if not (np.isfinite(ve).all() and np.isfinite(vn).all()):
        raise ValueError("Non-finite unmasked velocity predictions")
    lon, lat = np.meshgrid(longitude, latitude)
    distance = None
    supported = np.ones(lon.shape, dtype=bool)
    if args.mask_distance_km > 0:
        distance = nearest_distance_km(lon, lat, data)
        supported = distance <= args.mask_distance_km
    masked_e = np.where(supported, ve, np.nan)
    masked_n = np.where(supported, vn, np.nan)
    arrays = dict(longitude=longitude, latitude=latitude,
                        ve_raw_mm_yr=ve, vn_raw_mm_yr=vn, ve_mm_yr=masked_e, vn_mm_yr=masked_n,
                        support_mask=supported)
    if distance is not None:
        save_grid(gmt, output / "ve_eu_masked.nc", lon, lat, masked_e, args, env)
        save_grid(gmt, output / "vn_eu_masked.nc", lon, lat, masked_n, args, env)
        save_grid(gmt, output / "nearest_station_km.nc", lon, lat, distance, args, env, units="km")
        arrays["nearest_station_km"] = distance
    np.savez_compressed(output / "velocity_eu.npz", **arrays)
    # GMT -E reports direct model evaluations at the input sites, not samples
    # from the coarser 0.1-degree output grid. Inspect its documented layout.
    residual = np.loadtxt(output / "gmt_station_misfit.txt", comments="#", ndmin=2)
    if residual.shape[0] != len(data) or residual.shape[1] < 8:
        raise ValueError("Unexpected GMT -E residual layout; inspect gmt_station_misfit.txt")
    if not (np.allclose(residual[:, :2], data[:, :2], atol=1e-4) and
            np.allclose(residual[:, 2], data[:, 2], atol=1e-4) and
            np.allclose(residual[:, 5], data[:, 3], atol=1e-4)):
        raise ValueError("GMT residual columns do not match input; inspect raw residual file")
    # Columns: x y u_observed u_model u_residual v_observed v_model v_residual ...
    pred = residual[:, [3, 6]]
    error = data[:, 2:4] - pred
    block_mode = args.block_spacing is not None
    fit_name = "block_fit_eu.txt" if block_mode else "station_fit_eu.txt"
    np.savetxt(output / fit_name, np.column_stack((data, pred, error)), fmt="%.8f",
               header="lon lat VE_obs VN_obs sigmaE sigmaN VE_model VN_model VE_obs_minus_model VN_obs_minus_model; velocities: mm/yr")
    in_region = ((data[:, 0] >= args.region[0]) & (data[:, 0] <= args.region[1]) &
                 (data[:, 1] >= args.region[2]) & (data[:, 1] <= args.region[3]))
    statistics = {}
    for label, selection in ((("all_block_constraints" if block_mode else "all_input_stations"), np.ones(len(data), dtype=bool)),
                              (("output_region_blocks" if block_mode else "output_region_stations"), in_region)):
        err, sigma = error[selection], data[selection, 4:6]
        statistics[label] = {"count": int(selection.sum()),
                             "rms_E_mm_yr": float(np.sqrt(np.mean(err[:, 0]**2))),
                             "rms_N_mm_yr": float(np.sqrt(np.mean(err[:, 1]**2))),
                             "normalized_rms_E": float(np.sqrt(np.mean((err[:, 0]/sigma[:, 0])**2))),
                             "normalized_rms_N": float(np.sqrt(np.mean((err[:, 1]/sigma[:, 1])**2)))}
    statistics["grid_support"] = {"nearest_station_threshold_km": args.mask_distance_km,
                                   "fraction_visible": float(supported.mean()),
                                   "note": ("Mask disabled; complete Run 6.1-style grid" if distance is None else
                                            "Distance-based display mask, not a formal resolution or uncertainty estimate")}
    statistics["validation_status"] = ("Fit to block representative velocities only, not residuals at the original stations; no held-out validation"
                                       if block_mode else "In-sample fit only; no held-out validation or parameter optimization")
    if block_mode:
        statistics["input_kind"] = "blockmedian representatives; sigma is a representative scale"
    (output / "fit_summary.json").write_text(json.dumps(statistics, indent=2)+"\n")
    make_plot(output, args, lon, lat, masked_e, masked_n, data[in_region], error[in_region], statistics)
    return statistics


def make_plot(output, args, lon, lat, ve, vn, stations, error, statistics):
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(17, 11), layout="constrained")
    nu_label = "GMT default" if args.poisson is None else f"{args.poisson:g}"
    fudge_label = "GMT default" if args.fudge_km is None else f"Fd={args.fudge_km:g} km"
    block_label = f"blockmedian={args.block_spacing}, " if args.block_spacing else ""
    fig.suptitle(f"Run 3 | Run 6.1 gpsgridder method | Eurasia-fixed\n"
                  f"{block_label}"
                  f"retain={args.retain_percent:g}%, pixel grid={args.spacing}, "
                  f"nu={nu_label}, F={fudge_label}", fontsize=17)
    vmax = max(1, float(np.nanpercentile(np.abs(np.concatenate((ve.ravel(), vn.ravel()))), 99)))
    vnorm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    enorm = max(1, float(np.percentile(np.abs(error), 98)))
    residual_norm = TwoSlopeNorm(vmin=-enorm, vcenter=0, vmax=enorm)
    counts = statistics["output_region_blocks" if args.block_spacing else "output_region_stations"]
    for i, (ax, field, component) in enumerate(zip(axes[0], (ve, vn), ("East", "North"))):
        draw_context(ax, args.region, ROOT / "GNSS_analysis/data/basemap")
        im = ax.pcolormesh(lon, mercator_y(lat), field, cmap="RdBu_r", norm=vnorm, shading="auto", zorder=0.5, rasterized=True)
        ax.scatter(stations[:, 0], mercator_y(stations[:, 1]), s=1.5, c="#273c4b", alpha=0.35, zorder=2)
        mask_label = "full grid; no distance mask" if args.mask_distance_km == 0 else f"masked beyond {args.mask_distance_km:g} km"
        ax.set_title(f"{'ab'[i]}  {component} velocity | {mask_label}", loc="left", fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, extend="both", label="Velocity (mm/yr)")
    for i, (ax, component) in enumerate(zip(axes[1], ("East", "North"))):
        draw_context(ax, args.region, ROOT / "GNSS_analysis/data/basemap")
        im = ax.scatter(stations[:, 0], mercator_y(stations[:, 1]), c=error[:, i], s=10,
                         cmap="RdBu_r", norm=residual_norm, edgecolors="none", zorder=3)
        ax.set_title(f"{'cd'[i]}  Observed - fitted {component.lower()} | RMS={counts['rms_'+('E' if i == 0 else 'N')+'_mm_yr']:.3g} mm/yr", loc="left", fontsize=11)
        fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02, extend="both", label="In-sample residual (mm/yr)")
    fig.supxlabel("Run 6.1 solver settings. Training residuals are not independent validation; grid spacing is not strain resolution.", fontsize=10)
    fig.savefig(output / "run3_velocity_fit_eu.png", dpi=args.dpi, facecolor="white")
    fig.savefig(output / "run3_velocity_fit_eu.pdf", dpi=args.dpi, facecolor="white")
    plt.close(fig)

