from __future__ import annotations

import numpy as np


def temperature_level(T: np.ndarray, Tu: float, Tb: float, level: float = 0.5) -> float:
    return float(Tu + float(level) * (Tb - Tu))


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
        mask = np.isclose(r_sorted, rr, rtol=0.0, atol=max(1e-12, (bins[1] - bins[0]) * 0.5 if len(bins) > 1 else 1e-12))
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
