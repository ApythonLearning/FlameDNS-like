from __future__ import annotations

import numpy as np
import pandas as pd

from .ignition import centered_flame_coordinate
from .mesh import AxisymmetricMesh


def temperature_level(T: np.ndarray, Tu: float, Tb: float, level: float = 0.5) -> float:
    """Return the flame-front temperature level in K.

    The default ``level=0.5`` gives ``T_iso = (Tu + Tb) / 2``.
    """
    return float(Tu + float(level) * (Tb - Tu))


def midpoint_temperature(Tu: float, Tb: float) -> float:
    """Convenience wrapper for the requested ``(Tu + Tb) / 2`` isotherm."""
    return temperature_level(np.empty(0), Tu, Tb, level=0.5)


def temperature_field_from_profile(
    mesh: AxisymmetricMesh,
    profile: pd.DataFrame,
    flame_radius_m: float,
) -> np.ndarray:
    """Map a 1D temperature profile to a spherical front of radius ``R_f``.

    The 1D profile is centered at its maximum temperature gradient, then
    sampled at ``mesh.radius - R_f``. This mirrors the Stage 1 initialization
    and provides a temperature field from which the flame front is extracted.
    """
    x_rel = centered_flame_coordinate(profile)
    T = profile["T_K"].to_numpy(dtype=float)
    query = mesh.radius - float(flame_radius_m)
    return np.interp(query.ravel(), x_rel, T, left=T[0], right=T[-1]).reshape(mesh.radius.shape)


def extract_temperature_isoline(
    mesh: AxisymmetricMesh,
    T_K: np.ndarray,
    T_iso_K: float,
    time_s: float | None = None,
    flame_radius_m: float | None = None,
) -> pd.DataFrame:
    """Extract a 2D axisymmetric flame-front contour by marching squares.

    Output rows are ordered line-segment endpoints. Columns:
    ``segment_id, point_id, r_m, z_m, T_iso_K`` plus optional ``time_s`` and
    ``R_f_m``. Plot consecutive ``point_id=0,1`` rows with the same
    ``segment_id`` as one contour segment.
    """
    if T_K.shape != mesh.R.shape:
        raise ValueError(f"T_K shape {T_K.shape} does not match mesh shape {mesh.R.shape}.")

    records: list[dict[str, float | int]] = []
    segment_id = 0
    for i in range(len(mesh.r) - 1):
        for j in range(len(mesh.z) - 1):
            corners = [
                (mesh.r[i], mesh.z[j], T_K[i, j]),
                (mesh.r[i + 1], mesh.z[j], T_K[i + 1, j]),
                (mesh.r[i + 1], mesh.z[j + 1], T_K[i + 1, j + 1]),
                (mesh.r[i], mesh.z[j + 1], T_K[i, j + 1]),
            ]
            points = _cell_edge_crossings(corners, float(T_iso_K))
            if len(points) < 2:
                continue
            for start in range(0, len(points) - 1, 2):
                pair = points[start : start + 2]
                if len(pair) != 2:
                    continue
                for point_id, (r_m, z_m) in enumerate(pair):
                    row: dict[str, float | int] = {
                        "segment_id": segment_id,
                        "point_id": point_id,
                        "r_m": float(r_m),
                        "z_m": float(z_m),
                        "T_iso_K": float(T_iso_K),
                    }
                    if time_s is not None:
                        row["time_s"] = float(time_s)
                    if flame_radius_m is not None:
                        row["R_f_m"] = float(flame_radius_m)
                    records.append(row)
                segment_id += 1
    return pd.DataFrame.from_records(records)


def extract_front_history_from_profiles(
    mesh: AxisymmetricMesh,
    profile: pd.DataFrame,
    history: pd.DataFrame,
    Tu_K: float,
    Tb_K: float,
    sample_count: int = 6,
) -> pd.DataFrame:
    """Extract ``T=(Tu+Tb)/2`` fronts for selected times in a history table."""
    if history.empty:
        return pd.DataFrame(columns=["segment_id", "point_id", "r_m", "z_m", "T_iso_K", "time_s", "R_f_m"])

    sample_count = max(1, int(sample_count))
    indices = np.linspace(0, len(history) - 1, min(sample_count, len(history))).round().astype(int)
    T_iso = midpoint_temperature(Tu_K, Tb_K)
    frames = []
    next_segment_offset = 0
    for idx in np.unique(indices):
        row = history.iloc[int(idx)]
        R_f = float(row["R_f_m"])
        T = temperature_field_from_profile(mesh, profile, R_f)
        front = extract_temperature_isoline(mesh, T, T_iso, time_s=float(row["t_s"]), flame_radius_m=R_f)
        if front.empty:
            continue
        front["segment_id"] = front["segment_id"] + next_segment_offset
        next_segment_offset = int(front["segment_id"].max()) + 1
        frames.append(front)
    if not frames:
        return pd.DataFrame(columns=["segment_id", "point_id", "r_m", "z_m", "T_iso_K", "time_s", "R_f_m"])
    return pd.concat(frames, ignore_index=True)


def radius_from_temperature_isosurface(radius: np.ndarray, T: np.ndarray, T_iso: float) -> float:
    flat_r = radius.ravel()
    flat_T = T.ravel()
    idx = np.argsort(flat_r)
    r_sorted = flat_r[idx]
    T_sorted = flat_T[idx]
    bins = np.unique(r_sorted)
    if len(bins) > 2000:
        bins = np.linspace(r_sorted.min(), r_sorted.max(), 2000)
    radial_T = np.empty_like(bins)
    for i, rr in enumerate(bins):
        atol = max(1e-12, (bins[1] - bins[0]) * 0.5 if len(bins) > 1 else 1e-12)
        mask = np.isclose(r_sorted, rr, rtol=0.0, atol=atol)
        radial_T[i] = float(np.mean(T_sorted[mask])) if np.any(mask) else np.nan
    valid = np.isfinite(radial_T)
    bins = bins[valid]
    radial_T = radial_T[valid]
    diff = radial_T - T_iso
    crossing = np.where(np.signbit(diff[:-1]) != np.signbit(diff[1:]))[0]
    if len(crossing) == 0:
        return float(bins[int(np.argmin(np.abs(diff)))])
    i = int(crossing[0])
    return float(np.interp(T_iso, [radial_T[i], radial_T[i + 1]], [bins[i], bins[i + 1]]))


def radius_from_max_temperature_gradient(radius: np.ndarray, T: np.ndarray) -> float:
    flat_r = radius.ravel()
    flat_T = T.ravel()
    idx = np.argsort(flat_r)
    r = flat_r[idx]
    temp = flat_T[idx]
    unique_r, inverse = np.unique(r, return_inverse=True)
    radial_T = np.zeros_like(unique_r)
    counts = np.zeros_like(unique_r)
    np.add.at(radial_T, inverse, temp)
    np.add.at(counts, inverse, 1.0)
    radial_T /= np.maximum(counts, 1.0)
    grad = np.abs(np.gradient(radial_T, unique_r))
    return float(unique_r[int(np.argmax(grad))])


def _cell_edge_crossings(corners: list[tuple[float, float, float]], T_iso: float) -> list[tuple[float, float]]:
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    points: list[tuple[float, float]] = []
    for a, b in edges:
        r0, z0, T0 = corners[a]
        r1, z1, T1 = corners[b]
        d0 = float(T0) - T_iso
        d1 = float(T1) - T_iso
        if d0 == 0.0 and d1 == 0.0:
            continue
        if d0 == 0.0:
            point = (r0, z0)
        elif d1 == 0.0:
            point = (r1, z1)
        elif np.signbit(d0) == np.signbit(d1):
            continue
        else:
            weight = (T_iso - float(T0)) / (float(T1) - float(T0))
            point = (r0 + weight * (r1 - r0), z0 + weight * (z1 - z0))
        if not _contains_point(points, point):
            points.append(point)
    return points


def _contains_point(points: list[tuple[float, float]], point: tuple[float, float]) -> bool:
    return any(abs(point[0] - old[0]) < 1e-14 and abs(point[1] - old[1]) < 1e-14 for old in points)
