from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract flame radii from real PeleLMeX plotfiles.")
    parser.add_argument("plotfiles", nargs="+", help="PeleLMeX/AMReX plotfile directories, e.g. plt00020.")
    parser.add_argument("--summary", required=True, help="Cantera 1D summary JSON with T0_K and Tb_K.")
    parser.add_argument("--output-dir", default="postprocess_pelelmex")
    parser.add_argument("--temperature-field", default=None, help="Optional exact yt field name.")
    parser.add_argument("--iso-temperature", type=float, default=None, help="Optional temperature contour level in K.")
    parser.add_argument("--plot-last-field", action="store_true", help="Plot only the final plotfile temperature field.")
    parser.add_argument(
        "--last-field-margin-delta",
        type=float,
        default=2.5,
        help="Crop margin around the active final-time temperature field, in thermal flame thicknesses.",
    )
    parser.add_argument("--cfl", type=float, default=None, help="Optional run CFL for DNS warning output.")
    args = parser.parse_args()

    try:
        import yt
    except Exception as exc:
        raise SystemExit("yt is required to read PeleLMeX plotfiles. Install with `python -m pip install yt`.") from exc

    summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))
    Tu = float(summary["T0_K"])
    Tb = float(summary["Tb_K"])
    delta_T_m = float(summary["delta_T_m"])
    iso_temperature = args.iso_temperature if args.iso_temperature is not None else 0.5 * (Tu + Tb)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    contours_for_plot = []
    field_snapshots = []
    for plotfile in args.plotfiles:
        if not Path(plotfile).exists():
            raise SystemExit(f"Plotfile does not exist: {plotfile}")
        ds = yt.load(plotfile)
        if not rows:
            _print_dns_resolution_warning(ds, delta_T_m, args.cfl)
        field = _select_temperature_field(ds, args.temperature_field)
        r, z, temperature = _load_2d_temperature(ds, field)
        field_snapshots.append(
            (
                float(ds.current_time.to_value("s")) if hasattr(ds.current_time, "to_value") else float(ds.current_time),
                r,
                z,
                temperature,
            )
        )
        contour_segments = _contour_segments(r, z, temperature, iso_temperature)
        radii = _radii_from_contours(contour_segments)
        row = {
            "plotfile": str(plotfile),
            "time_s": float(ds.current_time.to_value("s")) if hasattr(ds.current_time, "to_value") else float(ds.current_time),
            "iso_temperature_K": iso_temperature,
            **radii,
        }
        rows.append(row)
        contours_for_plot.append((row["time_s"], contour_segments))

    history = pd.DataFrame(rows).sort_values("time_s")
    history.to_csv(output_dir / "flame_radii.csv", index=False)
    _plot_radius_history(history, output_dir / "flame_radii.png")
    _plot_contours(contours_for_plot, output_dir / "flame_contours.png")
    if args.plot_last_field:
        final_snapshot = max(field_snapshots, key=lambda item: item[0])
        _plot_last_temperature_field(
            final_snapshot,
            delta_T_m,
            iso_temperature,
            output_dir / "temperature_last.png",
            args.last_field_margin_delta,
        )
    print(history)


def _select_temperature_field(ds, requested: str | None):
    if requested:
        for field in ds.field_list + ds.derived_field_list:
            if str(field) == requested or field[-1] == requested:
                return field
        raise ValueError(f"Requested temperature field {requested!r} not found.")
    candidates = {"temp", "Temp", "temperature", "Temperature", "T"}
    for field in ds.field_list + ds.derived_field_list:
        if field[-1] in candidates:
            return field
    available = ", ".join(str(field) for field in ds.field_list[:30])
    raise ValueError(f"No temperature field found. First available fields: {available}")


def _print_dns_resolution_warning(ds, delta_T_m: float, cfl: float | None) -> None:
    dims = np.asarray(ds.domain_dimensions, dtype=float)
    width = np.asarray(ds.domain_width, dtype=float)
    active = dims > 1
    dx_min = float(np.min(width[active] / dims[active]))
    print(f"DNS check: dx_min={dx_min:.8e} m, delta_T={delta_T_m:.8e} m")
    if dx_min > delta_T_m / 10.0:
        print("WARNING: grid is under-resolved for DNS criterion dx <= delta_T/10.")
    elif dx_min > delta_T_m / 20.0:
        print("WARNING: grid passes dx <= delta_T/10 but not recommended dx <= delta_T/20.")
    if cfl is not None:
        print(f"DNS check: CFL={cfl:.4f}")
        if cfl >= 0.5:
            print("WARNING: CFL should be below 0.5.")


def _load_2d_temperature(ds, field):
    dims = np.array(ds.domain_dimensions, dtype=int)
    left = ds.domain_left_edge
    right = ds.domain_right_edge
    grid = ds.covering_grid(level=0, left_edge=left, dims=dims)
    values = np.asarray(grid[field], dtype=float)
    active_axes = [axis for axis, size in enumerate(values.shape) if size > 1]
    if len(active_axes) < 2:
        raise ValueError(f"Temperature field is not at least two-dimensional: shape={values.shape}")
    if len(active_axes) > 2:
        inactive = [axis for axis in range(values.ndim) if axis not in active_axes[:2]]
        for axis in sorted(inactive, reverse=True):
            values = np.take(values, indices=0, axis=axis)
    values = np.squeeze(values)

    axis0, axis1 = active_axes[:2]
    r = _cell_centers(float(left[axis0]), float(right[axis0]), int(dims[axis0]))
    z = _cell_centers(float(left[axis1]), float(right[axis1]), int(dims[axis1]))
    if values.shape != (len(r), len(z)):
        values = values.reshape((len(r), len(z)))
    return r, z, values


def _cell_centers(lo: float, hi: float, n: int) -> np.ndarray:
    dx = (hi - lo) / n
    return lo + dx * (np.arange(n) + 0.5)


def _contour_segments(r: np.ndarray, z: np.ndarray, temperature: np.ndarray, level: float) -> list[np.ndarray]:
    R, Z = np.meshgrid(r, z, indexing="ij")
    fig, ax = plt.subplots()
    contour = ax.contour(R, Z, temperature, levels=[level])
    segments: list[np.ndarray] = []
    if hasattr(contour, "allsegs"):
        for level_segments in contour.allsegs:
            for vertices in level_segments:
                if len(vertices) >= 2:
                    segments.append(np.asarray(vertices, dtype=float).copy())
    elif hasattr(contour, "collections"):
        for collection in contour.collections:
            for path in collection.get_paths():
                vertices = path.vertices
                if len(vertices) >= 2:
                    segments.append(vertices.copy())
    plt.close(fig)
    if not segments:
        t_min = float(np.nanmin(temperature))
        t_max = float(np.nanmax(temperature))
        raise ValueError(
            f"No T={level} K contour found in plotfile. "
            f"Temperature range is [{t_min}, {t_max}] K; "
            "choose a value inside this range with --iso-temperature."
        )
    return segments


def _radii_from_contours(segments: list[np.ndarray]) -> dict[str, float]:
    points = np.vstack(segments)
    radius = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    upper = radius[points[:, 1] >= 0.0]
    lower = radius[points[:, 1] <= 0.0]
    return {
        "radius_upper_m": float(np.max(upper)) if len(upper) else np.nan,
        "radius_lower_m": float(np.max(lower)) if len(lower) else np.nan,
        "radius_mean_m": float(np.mean(radius)),
        "radius_max_m": float(np.max(radius)),
        "contour_point_count": int(len(points)),
    }


def _plot_radius_history(history: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(history["time_s"], history["radius_upper_m"], label="upper")
    ax.plot(history["time_s"], history["radius_lower_m"], label="lower")
    ax.plot(history["time_s"], history["radius_mean_m"], label="mean")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("flame radius [m]")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_contours(contours_for_plot: list[tuple[float, list[np.ndarray]]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    for time_s, segments in contours_for_plot:
        for segment in segments:
            ax.plot(segment[:, 0], segment[:, 1], linewidth=1.0, label=f"{time_s:.4e} s")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("r [m]")
    ax.set_ylabel("z [m]")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _plot_last_temperature_field(
    snapshot: tuple[float, np.ndarray, np.ndarray, np.ndarray],
    delta_T_m: float,
    iso_temperature: float,
    path: Path,
    margin_delta: float,
) -> None:
    time_s, r, z, temperature = snapshot
    t_min = float(np.nanmin(temperature))
    t_max = float(np.nanmax(temperature))
    threshold = t_min + 0.05 * (t_max - t_min)
    active = temperature > threshold
    if np.any(active):
        active_r = np.where(np.any(active, axis=1))[0]
        active_z = np.where(np.any(active, axis=0))[0]
        margin_m = margin_delta * delta_T_m
        r_lo = max(float(r[0]), float(r[active_r[0]]) - margin_m)
        r_hi = min(float(r[-1]), float(r[active_r[-1]]) + margin_m)
        z_lo = max(float(z[0]), float(z[active_z[0]]) - margin_m)
        z_hi = min(float(z[-1]), float(z[active_z[-1]]) + margin_m)
    else:
        r_lo, r_hi = float(r[0]), float(r[-1])
        z_lo, z_hi = float(z[0]), float(z[-1])

    r_mask = (r >= r_lo) & (r <= r_hi)
    z_mask = (z >= z_lo) & (z <= z_hi)
    r_plot = r[r_mask] / delta_T_m
    z_plot = z[z_mask] / delta_T_m
    temperature_plot = temperature[np.ix_(r_mask, z_mask)].T

    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    mesh = ax.pcolormesh(
        r_plot,
        z_plot,
        temperature_plot,
        shading="auto",
        cmap="coolwarm",
        vmin=t_min,
        vmax=t_max,
    )
    if t_min < iso_temperature < t_max:
        ax.contour(
            r_plot,
            z_plot,
            temperature_plot,
            levels=[iso_temperature],
            colors="white",
            linewidths=0.8,
        )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$x/\delta_T$")
    ax.set_ylabel(r"$y/\delta_T$")
    ax.text(
        0.98,
        0.94,
        rf"$t = {time_s:.3e}\ \mathrm{{s}}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11,
        fontweight="bold",
    )
    fig.colorbar(mesh, ax=ax, label="T [K]", fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(path, dpi=240)
    plt.close(fig)


if __name__ == "__main__":
    main()
