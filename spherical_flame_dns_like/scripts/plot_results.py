from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.flame_front import extract_front_history_from_profiles
from src.mesh import create_axisymmetric_mesh
from src.utils import ensure_dir, load_yaml


def save_figure(fig: plt.Figure, out_base: Path) -> None:
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=200)
    fig.savefig(out_base.with_suffix(".svg"))
    plt.close(fig)


def plot_front_history(front: pd.DataFrame, case_name: str, figures_dir: Path) -> None:
    if front.empty:
        return

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    times = sorted(front["time_s"].unique())
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(min(times), max(times) if len(times) > 1 else min(times) + 1.0)

    for time_s in times:
        subset = front[front["time_s"] == time_s]
        segments = []
        mirrored_segments = []
        for _, segment in subset.groupby("segment_id"):
            ordered = segment.sort_values("point_id")
            if len(ordered) != 2:
                continue
            r = ordered["r_m"].to_numpy() * 1000.0
            z = ordered["z_m"].to_numpy() * 1000.0
            segments.append(list(zip(r, z)))
            mirrored_segments.append(list(zip(-r, z)))
        color = cmap(norm(time_s))
        ax.add_collection(LineCollection(segments, colors=[color], linewidths=1.4))
        ax.add_collection(LineCollection(mirrored_segments, colors=[color], linewidths=0.9, alpha=0.45))

    scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(scalar_map, ax=ax, pad=0.02)
    cbar.set_label("t [s]")
    ax.axvline(0.0, color="0.25", linewidth=0.8)
    max_r = float(front["r_m"].max() * 1000.0)
    max_z = float(front["z_m"].abs().max() * 1000.0)
    lim = max(max_r, max_z) * 1.08
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("r [mm], mirrored for visualization")
    ax.set_ylabel("z [mm]")
    ax.set_title(f"{case_name}: T=(Tu+Tb)/2 flame-front evolution")
    ax.grid(True, alpha=0.25)
    save_figure(fig, figures_dir / f"flame_front_evolution_{case_name}")


def plot_sb_kappa(summary: pd.DataFrame, results_dir: Path, figures_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    colors = plt.cm.plasma(np.linspace(0.08, 0.92, len(summary)))
    markers = ["*", "p", "D", "d", "^", "o", "s"]

    for idx, (_, row) in enumerate(summary.sort_values("h2_volume_fraction").iterrows()):
        history = pd.read_csv(results_dir / row["case"] / "flame_history.csv")
        kappa = history["kappa_1_s"].to_numpy()
        Sb = history["S_b_m_s"].to_numpy()
        color = colors[idx]
        marker = markers[idx % len(markers)]
        label = f"{100 * row['h2_volume_fraction']:.0f}% H2"
        ax.scatter(
            kappa,
            Sb,
            s=28,
            marker=marker,
            facecolors="none" if marker not in ["*", "."] else color,
            edgecolors=color,
            linewidths=1.2,
            alpha=0.9,
            label=label,
        )
        if len(kappa) >= 3:
            slope, intercept = np.polyfit(kappa, Sb, 1)
            x_fit = np.linspace(float(kappa.min()), float(kappa.max()), 120)
            ax.plot(x_fit, intercept + slope * x_fit, linestyle=":", linewidth=1.8, color=color)

    ax.set_xlabel("stretch rate kappa [1/s]")
    ax.set_ylabel("burning velocity S_b [m/s]")
    ax.legend(ncol=1, fontsize=8, frameon=False)
    ax.grid(True, alpha=0.28)
    ax.set_title("S_b-kappa relation with linear fits")
    save_figure(fig, figures_dir / "Sb_kappa")


def main() -> None:
    config = load_yaml(ROOT / "config" / "base.yaml")
    results_dir = ROOT / config["project"]["output_dir"]
    figures_dir = ensure_dir(results_dir / "figures")
    summary_path = results_dir / "sweep_summary.csv"
    if not summary_path.exists():
        raise SystemExit("Missing results/sweep_summary.csv. Run python scripts/run_sweep.py first.")
    summary = pd.read_csv(summary_path)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for _, row in summary.iterrows():
        history = pd.read_csv(results_dir / row["case"] / "flame_history.csv")
        ax.plot(history["t_s"] * 1000.0, history["R_f_m"] * 1000.0, label=f"{100 * row['h2_volume_fraction']:.0f}% H2")
    ax.set_xlabel("t [ms]")
    ax.set_ylabel("R_f [mm]")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Rf_t")

    plot_sb_kappa(summary, results_dir, figures_dir)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for _, row in summary.iterrows():
        history = pd.read_csv(results_dir / row["case"] / "flame_history.csv")
        ax.plot(history["R_f_m"] * 1000.0, history["S_b_m_s"], label=f"{100 * row['h2_volume_fraction']:.0f}% H2")
    ax.set_xlabel("R_f [mm]")
    ax.set_ylabel("S_b [m/s]")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Sb_Rf")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(summary["h2_volume_fraction"] * 100.0, summary["Ma"], "o-")
    ax.set_xlabel("H2 volume fraction [%]")
    ax.set_ylabel("Markstein number Ma = L_b / delta_f [-]")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Ma_h2")

    front_case = summary.sort_values("h2_volume_fraction").iloc[-1]
    case_dir = results_dir / front_case["case"]
    profile = pd.read_csv(case_dir / "free_flame_profile.csv")
    history = pd.read_csv(case_dir / "flame_history.csv")
    mesh = create_axisymmetric_mesh(config["mesh"])
    sample_count = int(config.get("postprocess", {}).get("front_sample_count", 7))
    front = extract_front_history_from_profiles(
        mesh=mesh,
        profile=profile,
        history=history,
        Tu_K=float(front_case["T_unburned_K"]),
        Tb_K=float(front_case["T_burned_K"]),
        sample_count=sample_count,
    )
    front.to_csv(case_dir / "flame_front_history.csv", index=False)
    plot_front_history(front, str(front_case["case"]), figures_dir)

    print(f"Saved figures to {figures_dir}")


if __name__ == "__main__":
    main()
