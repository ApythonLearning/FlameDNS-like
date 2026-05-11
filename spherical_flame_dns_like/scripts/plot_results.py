from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.utils import ensure_dir, load_yaml


def save_figure(fig: plt.Figure, out_base: Path) -> None:
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=200)
    fig.savefig(out_base.with_suffix(".svg"))
    plt.close(fig)


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

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for _, row in summary.iterrows():
        history = pd.read_csv(results_dir / row["case"] / "flame_history.csv")
        ax.plot(history["kappa_1_s"], history["S_b_m_s"], ".", ms=3, label=f"{100 * row['h2_volume_fraction']:.0f}% H2")
    ax.set_xlabel("kappa [1/s]")
    ax.set_ylabel("S_b [m/s]")
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Sb_kappa")

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

    print(f"Saved figures to {figures_dir}")


if __name__ == "__main__":
    main()
