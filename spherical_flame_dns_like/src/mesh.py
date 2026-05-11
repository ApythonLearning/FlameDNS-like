from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class AxisymmetricMesh:
    r: np.ndarray
    z: np.ndarray
    R: np.ndarray
    Z: np.ndarray
    radius: np.ndarray
    dr: float
    dz: float


def create_axisymmetric_mesh(mesh_config: dict) -> AxisymmetricMesh:
    r = np.linspace(0.0, float(mesh_config["r_max_m"]), int(mesh_config["nr"]))
    z = np.linspace(float(mesh_config["z_min_m"]), float(mesh_config["z_max_m"]), int(mesh_config["nz"]))
    R, Z = np.meshgrid(r, z, indexing="ij")
    radius = np.sqrt(R**2 + Z**2)
    dr = float(r[1] - r[0]) if len(r) > 1 else 0.0
    dz = float(z[1] - z[0]) if len(z) > 1 else 0.0
    return AxisymmetricMesh(r=r, z=z, R=R, Z=Z, radius=radius, dr=dr, dz=dz)


def points_per_flame_thickness(mesh: AxisymmetricMesh, delta_f_m: float) -> float:
    dx = min(mesh.dr, mesh.dz)
    if dx <= 0.0:
        return 0.0
    return float(delta_f_m / dx)
