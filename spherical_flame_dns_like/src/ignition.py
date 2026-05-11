from __future__ import annotations

import numpy as np
import pandas as pd

from .mesh import AxisymmetricMesh


def centered_flame_coordinate(profile: pd.DataFrame) -> np.ndarray:
    x = profile["x_m"].to_numpy()
    T = profile["T_K"].to_numpy()
    grad_idx = int(np.argmax(np.abs(np.gradient(T, x))))
    return x - x[grad_idx]


def map_profile_to_axisymmetric_mesh(
    mesh: AxisymmetricMesh,
    profile: pd.DataFrame,
    initial_radius_m: float,
) -> dict[str, np.ndarray]:
    x_rel = centered_flame_coordinate(profile)
    query = mesh.radius - float(initial_radius_m)
    fields: dict[str, np.ndarray] = {}
    for column in profile.columns:
        if column == "x_m":
            continue
        values = profile[column].to_numpy()
        fields[column] = np.interp(query.ravel(), x_rel, values, left=values[0], right=values[-1]).reshape(mesh.radius.shape)
    return fields


def apply_hot_kernel(fields: dict[str, np.ndarray], mesh: AxisymmetricMesh, ignition_config: dict) -> dict[str, np.ndarray]:
    radius = float(ignition_config["hot_kernel_radius_m"])
    hot_T = float(ignition_config["hot_kernel_temperature_K"])
    mask = mesh.radius <= radius
    out = {key: value.copy() for key, value in fields.items()}
    if "T_K" in out:
        out["T_K"][mask] = np.maximum(out["T_K"][mask], hot_T)
    return out


def field_to_dataframe(mesh: AxisymmetricMesh, fields: dict[str, np.ndarray]) -> pd.DataFrame:
    data = {
        "r_m": mesh.R.ravel(),
        "z_m": mesh.Z.ravel(),
        "radius_m": mesh.radius.ravel(),
    }
    for key, value in fields.items():
        data[key] = value.ravel()
    return pd.DataFrame(data)
