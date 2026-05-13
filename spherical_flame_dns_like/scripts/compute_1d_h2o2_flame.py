from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.chemistry import equivalence_ratio_h2_o2
from src.utils import ensure_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute a 1D premixed H2/O2 free flame with Cantera.")
    parser.add_argument("--h2-volume-fraction", type=float, default=0.10)
    parser.add_argument("--phi", type=float, default=None, help="Overrides --h2-volume-fraction if provided.")
    parser.add_argument("--T0", type=float, default=298.0)
    parser.add_argument("--P0", type=float, default=101325.0)
    parser.add_argument("--mechanism", default="h2o2.yaml")
    parser.add_argument("--width", type=float, default=0.04)
    parser.add_argument("--transport", default="mixture-averaged")
    parser.add_argument("--output-dir", default=str(ROOT / "profiles" / "h2o2_phi_auto"))
    parser.add_argument("--loglevel", type=int, default=1)
    parser.add_argument("--skip-hdf5", action="store_true", help="Write CSV/JSON only. Use only for dependency debugging.")
    args = parser.parse_args()

    import cantera as ct

    composition, phi, h2_fraction = _composition(args.phi, args.h2_volume_fraction)
    output_dir = ensure_dir(args.output_dir)

    gas = ct.Solution(args.mechanism)
    gas.TPX = args.T0, args.P0, composition
    rho_u = float(gas.density)

    flame = ct.FreeFlame(gas, width=args.width)
    flame.transport_model = args.transport
    flame.set_refine_criteria(ratio=3.0, slope=0.04, curve=0.08)
    flame.solve(loglevel=args.loglevel, auto=True)

    x = np.asarray(flame.grid, dtype=float)
    T = np.asarray(flame.T, dtype=float)
    velocity = np.asarray(flame.velocity, dtype=float)
    density = np.asarray(flame.density, dtype=float)
    data: dict[str, np.ndarray] = {
        "x_m": x,
        "T_K": T,
        "u_m_s": velocity,
        "rho_kg_m3": density,
    }
    for species in gas.species_names:
        data[f"Y_{species}"] = np.asarray(flame.Y[gas.species_index(species), :], dtype=float)
        data[f"X_{species}"] = np.asarray(flame.X[gas.species_index(species), :], dtype=float)

    profile = pd.DataFrame(data)
    delta_T_m = _thermal_thickness(x, T)
    summary = {
        "mechanism": args.mechanism,
        "transport": args.transport,
        "phi": phi,
        "h2_volume_fraction": h2_fraction,
        "T0_K": float(args.T0),
        "P0_Pa": float(args.P0),
        "S_L_m_s": abs(float(velocity[0])),
        "delta_T_m": delta_T_m,
        "Tb_K": float(T[-1]),
        "rho_u_kg_m3": rho_u,
        "rho_b_kg_m3": float(density[-1]),
        "grid_points": int(len(x)),
    }

    profile.to_csv(output_dir / "h2o2_1d_profile.csv", index=False)
    pd.DataFrame([summary]).to_csv(output_dir / "h2o2_1d_summary.csv", index=False)
    (output_dir / "h2o2_1d_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if not args.skip_hdf5:
        _write_hdf5(output_dir / "h2o2_1d_profile.h5", profile, summary)
    print(json.dumps(summary, indent=2))


def _composition(phi: float | None, h2_volume_fraction: float) -> tuple[dict[str, float], float, float]:
    if phi is not None:
        if phi <= 0.0:
            raise ValueError("phi must be positive.")
        h2 = 2.0 * phi / (1.0 + 2.0 * phi)
        return {"H2": h2, "O2": 1.0 - h2}, float(phi), float(h2)
    if not 0.0 < h2_volume_fraction < 1.0:
        raise ValueError("h2-volume-fraction must be between 0 and 1.")
    return (
        {"H2": h2_volume_fraction, "O2": 1.0 - h2_volume_fraction},
        equivalence_ratio_h2_o2(h2_volume_fraction),
        float(h2_volume_fraction),
    )


def _thermal_thickness(x: np.ndarray, T: np.ndarray) -> float:
    dTdx = np.gradient(T, x)
    max_grad = float(np.max(np.abs(dTdx)))
    if max_grad <= 0.0:
        raise ValueError("Cannot compute thermal flame thickness from zero temperature gradient.")
    return float((T.max() - T.min()) / max_grad)


def _write_hdf5(path: Path, profile: pd.DataFrame, summary: dict[str, float | int | str]) -> None:
    try:
        import h5py
    except Exception as exc:
        raise SystemExit(
            "h5py is required for HDF5 profile output. Install dependencies with "
            "`python -m pip install -r requirements.txt`."
        ) from exc
    with h5py.File(path, "w") as handle:
        group = handle.create_group("profile")
        for column in profile.columns:
            group.create_dataset(column, data=profile[column].to_numpy())
        for key, value in summary.items():
            handle.attrs[key] = value


if __name__ == "__main__":
    main()
