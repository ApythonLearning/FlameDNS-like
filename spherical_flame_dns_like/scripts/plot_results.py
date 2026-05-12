from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.flame_front import extract_front_history_from_profiles
from src.mesh import create_axisymmetric_mesh
from src.utils import ensure_dir, load_yaml


def configure_reference_style() -> None:
    """Use a journal-like Matplotlib style close to the reference figures."""

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
            "axes.linewidth": 0.9,
            "axes.labelsize": 18,
            "axes.titlesize": 18,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 15,
            "legend.frameon": False,
            "legend.handlelength": 2.2,
            "legend.borderaxespad": 0.3,
            "savefig.bbox": "tight",
        }
    )


def apply_axis_style(ax: plt.Axes) -> None:
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.minorticks_on()
    ax.grid(False)


def apply_figure_style(fig: plt.Figure) -> None:
    for ax in fig.axes:
        apply_axis_style(ax)


def add_legend(ax: plt.Axes, **kwargs) -> None:
    defaults = {"frameon": False}
    defaults.update(kwargs)
    ax.legend(**defaults)


def add_panel_label(ax: plt.Axes, label: str, x: float = 0.04, y: float = 0.88) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=15,
        fontstyle="italic",
        va="top",
        ha="left",
    )


def add_le_eff_annotation_row(
    ax: plt.Axes,
    x_values: np.ndarray,
    le_eff_values: np.ndarray,
    y: float = 0.96,
) -> None:
    """Place Le_eff numeric values inside the plotting area instead of x tick labels."""

    transform = blended_transform_factory(ax.transData, ax.transAxes)
    ax.margins(y=0.22)
    for x_value, le_eff in zip(x_values, le_eff_values):
        ax.text(
            x_value,
            y,
            f"{le_eff:.3f}",
            transform=transform,
            fontsize=9,
            ha="center",
            va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.0},
        )


def save_figure(fig: plt.Figure, out_base: Path) -> None:
    apply_figure_style(fig)
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
    cbar.set_label(r"$t$ (s)")
    ax.axvline(0.0, color="0.25", linewidth=0.8)
    max_r = float(front["r_m"].max() * 1000.0)
    max_z = float(front["z_m"].abs().max() * 1000.0)
    lim = max(max_r, max_z) * 1.08
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$r$ (mm), mirrored")
    ax.set_ylabel(r"$z$ (mm)")
    ax.set_title(rf"{case_name}: $T=(T_u+T_b)/2$ flame front")
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
        label = rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$"
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

    ax.set_xlabel(r"$\kappa$ (s$^{-1}$)")
    ax.set_ylabel(r"$S_b$ (m s$^{-1}$)")
    add_legend(ax, ncol=1)
    ax.grid(True, alpha=0.28)
    ax.set_title(r"$S_b$-$\kappa$ relation")
    save_figure(fig, figures_dir / "Sb_kappa")


def plot_cellular_summary(summary: pd.DataFrame, figures_dir: Path) -> None:
    modes = [
        ("unstable", "RT-Unstable", "tab:red", "o"),
        ("neutral", "RT-neutral", "0.25", "s"),
        ("stable", "RT-Stable", "tab:blue", "^"),
    ]
    required = {"cellular_Le_eff", "cellular_rt_unstable_max_growth_1_s"}
    if not required.issubset(summary.columns):
        return

    ordered = summary.sort_values("h2_volume_fraction")
    x = ordered["h2_volume_fraction"].to_numpy() * 100.0
    xtick_labels = [rf"{100.0 * row.h2_volume_fraction:.0f}%" for row in ordered.itertuples()]
    le_eff_values = ordered["cellular_Le_eff"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for mode, label, color, marker in modes:
        column = f"cellular_rt_{mode}_max_growth_1_s"
        if column in ordered.columns:
            ax.plot(x, ordered[column], marker=marker, color=color, linewidth=1.7, label=label)
    ax.axhline(0.0, color="0.25", linewidth=0.9)
    ax.set_xlabel(r"$\mathrm{H}_2$ volume fraction (%)")
    ax.set_ylabel(r"$\omega_{\max}$ (s$^{-1}$)")
    ax.set_title(r"DL/TD/RT growth")
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels)
    add_le_eff_annotation_row(ax, x, le_eff_values)
    add_legend(ax, loc="best")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "cellular_growth_rt_modes")

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    for mode, label, color, marker in modes:
        column = f"cellular_rt_{mode}_lambda_over_delta_at_max"
        if column in ordered.columns:
            ax.plot(x, ordered[column], marker=marker, color=color, linewidth=1.7, label=label)
    ax.set_xlabel(r"$\mathrm{H}_2$ volume fraction (%)")
    ax.set_ylabel(r"$\lambda_{\max}/\delta_f$")
    ax.set_title(r"Most amplified wavelength")
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels)
    add_le_eff_annotation_row(ax, x, le_eff_values)
    add_legend(ax, loc="best")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "cellular_lambda_rt_modes")


def plot_dispersion_relations(summary: pd.DataFrame, results_dir: Path, figures_dir: Path) -> None:
    if "cellular_enabled" not in summary.columns:
        return

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(summary)))
    linestyles = {"unstable": "-", "neutral": "--", "stable": ":"}
    mode_labels = {"unstable": "RT-Unstable", "neutral": "RT-neutral", "stable": "RT-Stable"}
    for idx, (_, row) in enumerate(summary.sort_values("h2_volume_fraction").iterrows()):
        path = results_dir / row["case"] / "cellular_instability.csv"
        if not path.exists():
            continue
        dispersion = pd.read_csv(path)
        if dispersion.empty:
            continue
        final = dispersion[dispersion["t_s"] == dispersion["t_s"].max()]
        le_eff = float(final["Le_eff"].iloc[0]) if "Le_eff" in final else np.nan
        for mode, mode_data in final.groupby("rt_mode"):
            ax.plot(
                mode_data["k_delta"],
                mode_data["omega_dimensionless"],
                color=colors[idx],
                linestyle=linestyles.get(str(mode), "-"),
                linewidth=1.45,
                label=rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$, "
                rf"$Le_{{\rm eff}}={le_eff:.3f}$, {mode_labels.get(str(mode), str(mode))}",
            )
    ax.axhline(0.0, color="0.25", linewidth=0.9)
    ax.set_xlabel(r"$k\delta_f$")
    ax.set_ylabel(r"$\omega\delta_f/S_L$")
    ax.set_title(r"Law-model dispersion")
    add_legend(ax, ncol=2, fontsize=8.5)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "cellular_dispersion_rt_modes")


def plot_dispersion_grid_rt_modes(summary: pd.DataFrame, results_dir: Path, figures_dir: Path) -> None:
    if "cellular_enabled" not in summary.columns:
        return

    ordered = summary.sort_values("h2_volume_fraction").head(8).reset_index(drop=True)
    if ordered.empty:
        return

    mode_styles = {
        "unstable": ("RT-Unstable", "tab:red", "-"),
        "neutral": ("RT-neutral", "0.25", "--"),
        "stable": ("RT-Stable", "tab:blue", ":"),
    }
    fig, axes = plt.subplots(2, 4, figsize=(16.0, 8.2), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for panel_idx, (ax, (_, row)) in enumerate(zip(axes_flat, ordered.iterrows())):
        path = results_dir / row["case"] / "cellular_instability.csv"
        if not path.exists():
            ax.set_axis_off()
            continue
        dispersion = pd.read_csv(path)
        if dispersion.empty:
            ax.set_axis_off()
            continue
        final = dispersion[dispersion["t_s"] == dispersion["t_s"].max()]
        le_eff = float(final["Le_eff"].iloc[0]) if "Le_eff" in final else np.nan
        for mode, mode_data in final.groupby("rt_mode"):
            label, color, linestyle = mode_styles.get(str(mode), (str(mode), "0.4", "-"))
            ax.plot(
                mode_data["k_delta"],
                mode_data["omega_dimensionless"],
                color=color,
                linestyle=linestyle,
                linewidth=1.8,
                label=label,
            )
        ax.axhline(0.0, color="0.55", linewidth=0.8)
        add_panel_label(ax, f"({chr(ord('a') + panel_idx)})")
        ax.set_title(
            rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$, $Le_{{\rm eff}}={le_eff:.3f}$",
            fontsize=13,
            pad=8,
        )
        ax.grid(True, alpha=0.25)

    for ax in axes_flat[len(ordered) :]:
        ax.set_axis_off()
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$k\delta_f$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\omega\delta_f/S_L$")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02), fontsize=11)
    fig.suptitle(r"DL/TD/RT dispersion comparison", y=1.075, fontsize=15)
    save_figure(fig, figures_dir / "cellular_dispersion_2x4_rt_modes")


def plot_nonlinear_grid_rt_modes(summary: pd.DataFrame, results_dir: Path, figures_dir: Path) -> None:
    if "nonlinear_cellular_enabled" not in summary.columns:
        return

    ordered = summary.sort_values("h2_volume_fraction").head(8).reset_index(drop=True)
    if ordered.empty:
        return

    mode_styles = {
        "unstable": ("RT-Unstable", "tab:red", "-"),
        "neutral": ("RT-neutral", "0.25", "--"),
        "stable": ("RT-Stable", "tab:blue", ":"),
    }
    fig, axes = plt.subplots(2, 4, figsize=(16.0, 8.2), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for panel_idx, (ax, (_, row)) in enumerate(zip(axes_flat, ordered.iterrows())):
        path = results_dir / row["case"] / "nonlinear_cellular.csv"
        if not path.exists():
            ax.set_axis_off()
            continue
        nonlinear = pd.read_csv(path)
        if nonlinear.empty:
            ax.set_axis_off()
            continue
        le_eff = float(nonlinear["Le_eff"].iloc[0]) if "Le_eff" in nonlinear else np.nan
        threshold = float(nonlinear["nonlinear_threshold_over_delta"].iloc[0])
        for mode, mode_data in nonlinear.groupby("rt_mode"):
            label, color, linestyle = mode_styles.get(str(mode), (str(mode), "0.4", "-"))
            ax.plot(
                mode_data["t_s"] * 1000.0,
                mode_data["amplitude_over_delta"],
                color=color,
                linestyle=linestyle,
                linewidth=1.8,
                label=label,
            )
        ax.axhline(threshold, color="0.55", linewidth=0.8)
        add_panel_label(ax, f"({chr(ord('a') + panel_idx)})")
        ax.set_title(
            rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$, $Le_{{\rm eff}}={le_eff:.3f}$",
            fontsize=13,
            pad=8,
        )
        ax.grid(True, alpha=0.25)

    for ax in axes_flat[len(ordered) :]:
        ax.set_axis_off()
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$t$ (ms)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$A/\delta_f$")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.02), fontsize=11)
    fig.suptitle(r"Weakly nonlinear cellular-amplitude evolution", y=1.075, fontsize=15)
    save_figure(fig, figures_dir / "nonlinear_cellular_amplitude_2x4_rt_modes")


def plot_nonlinear_summary(summary: pd.DataFrame, figures_dir: Path) -> None:
    modes = [
        ("unstable", "RT-Unstable", "tab:red", "o"),
        ("neutral", "RT-neutral", "0.25", "s"),
        ("stable", "RT-Stable", "tab:blue", "^"),
    ]
    if "nonlinear_cellular_enabled" not in summary.columns:
        return

    ordered = summary.sort_values("h2_volume_fraction")
    x = ordered["h2_volume_fraction"].to_numpy() * 100.0
    xtick_labels = [rf"{100.0 * row.h2_volume_fraction:.0f}%" for row in ordered.itertuples()]
    le_eff_values = ordered["cellular_Le_eff"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for mode, label, color, marker in modes:
        column = f"nonlinear_rt_{mode}_onset_time_s"
        if column in ordered.columns:
            y = ordered[column].to_numpy(dtype=float) * 1000.0
            ax.plot(x, y, marker=marker, color=color, linewidth=1.7, label=label)
    ax.set_xlabel(r"$\mathrm{H}_2$ volume fraction (%)")
    ax.set_ylabel(r"$t_{\rm onset}$ (ms)")
    ax.set_title(r"Weakly nonlinear cellular onset")
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels)
    add_le_eff_annotation_row(ax, x, le_eff_values)
    add_legend(ax, loc="best")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "nonlinear_cellular_onset_rt_modes")

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for mode, label, color, marker in modes:
        column = f"nonlinear_rt_{mode}_cellular_speed_factor_final"
        if column in ordered.columns:
            ax.plot(x, ordered[column], marker=marker, color=color, linewidth=1.7, label=label)
    ax.set_xlabel(r"$\mathrm{H}_2$ volume fraction (%)")
    ax.set_ylabel(r"$S_{b,\rm cell}/S_b$")
    ax.set_title(r"Nonlinear cellular speed enhancement")
    ax.set_xticks(x)
    ax.set_xticklabels(xtick_labels)
    add_le_eff_annotation_row(ax, x, le_eff_values)
    add_legend(ax, loc="best")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "nonlinear_cellular_speed_factor_rt_modes")


def _time_sample_rows(frame: pd.DataFrame, sample_count: int) -> pd.DataFrame:
    if frame.empty:
        return frame
    ordered = frame.sort_values("t_s").reset_index(drop=True)
    indices = np.unique(np.linspace(0, len(ordered) - 1, max(sample_count, 2)).round().astype(int))
    return ordered.iloc[indices].reset_index(drop=True)


def plot_nonlinear_front_locations(summary: pd.DataFrame, results_dir: Path, figures_dir: Path, config: dict) -> None:
    if "nonlinear_cellular_enabled" not in summary.columns:
        return

    front_config = config.get("cellular_instability", {}).get("nonlinear", {}).get("front_visualization", {})
    case_name = str(front_config.get("case", ""))
    if case_name:
        selected = summary[summary["case"] == case_name]
        if selected.empty:
            return
        row = selected.iloc[0]
    else:
        row = summary.sort_values("h2_volume_fraction").iloc[-1]

    path = results_dir / row["case"] / "nonlinear_cellular.csv"
    if not path.exists():
        return
    nonlinear = pd.read_csv(path)
    if nonlinear.empty:
        return

    sample_count = int(front_config.get("sample_count", 6))
    y_max = float(front_config.get("y_over_delta_max", 60.0))
    amplitude_scale = float(front_config.get("amplitude_visual_scale", 160.0))
    base_shift = float(front_config.get("x_over_delta_shift", 18.0))
    y = np.linspace(0.0, y_max, 500)

    mode_order = [
        ("stable", "RT-Stable", "(a)"),
        ("neutral", "RT-neutral", "(b)"),
        ("unstable", "RT-Unstable", "(c)"),
    ]
    colors = plt.cm.tab10(np.linspace(0.0, 0.9, sample_count))
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.4), sharex=True, sharey=True)

    for ax, (mode, mode_label, panel_label) in zip(axes, mode_order):
        mode_data = nonlinear[nonlinear["rt_mode"] == mode].copy()
        if mode_data.empty:
            ax.set_axis_off()
            continue
        gravity = float(mode_data["gravity_effective_m_s2"].iloc[0])
        samples = _time_sample_rows(mode_data, sample_count)
        for color, (_, sample) in zip(colors, samples.iterrows()):
            lambda_over_delta = max(float(sample["dominant_lambda_over_delta"]), 1.0e-12)
            phase = 2.0 * np.pi * y / lambda_over_delta
            amplitude = amplitude_scale * float(sample["amplitude_over_delta"])
            x_center = base_shift + float(sample["R_f_m"] / sample["dominant_lambda_m"] * lambda_over_delta)

            if mode == "stable":
                # RT-stable gravity damps cellular deformation; keep only a weak,
                # broad displacement so this panel remains visibly distinct from
                # the RT-neutral sinusoidal cellular front.
                broad_shape = 2.0 * y / max(y_max, 1.0) - 1.0
                weak_cell = 0.10 * np.sin(phase + 0.35)
                perturbation = 0.08 * amplitude * (0.30 * broad_shape + weak_cell)
            elif mode == "neutral":
                perturbation = amplitude * (0.85 * np.sin(phase) + 0.15 * np.sin(2.0 * phase))
            else:
                center_envelope = np.exp(-((y - 0.5 * y_max) / (0.19 * y_max)) ** 2)
                gravity_sign = -1.0 if gravity < 0.0 else 1.0
                perturbation = gravity_sign * amplitude * (1.0 + 3.8 * center_envelope) * np.sin(phase)

            ax.plot(x_center + perturbation, y, color=color, linewidth=1.45)

        add_panel_label(ax, panel_label, x=0.03, y=0.86)
        ax.text(0.98, 0.82, mode_label, transform=ax.transAxes, fontsize=12, ha="right")
        # Sign convention follows cellular_instability.effective_gravity:
        # g > 0 is RT-stable and points left in this front-location plot,
        # g < 0 is RT-unstable and points right, matching the reference layout.
        if gravity > 0.0:
            ax.annotate("", xy=(10.0, 30.0), xytext=(28.0, 30.0), arrowprops={"arrowstyle": "simple", "color": "0.15"})
        elif gravity < 0.0:
            ax.annotate("", xy=(28.0, 30.0), xytext=(10.0, 30.0), arrowprops={"arrowstyle": "simple", "color": "0.15"})
        ax.set_ylabel(r"$y/\delta_f$")
        ax.set_ylim(0.0, y_max)
        ax.grid(False)

    axes[-1].set_xlabel(r"$x/\delta_f$")
    x_max = float(front_config.get("x_over_delta_max", 205.0))
    for ax in axes:
        ax.set_xlim(0.0, x_max)
    le_eff = float(nonlinear["Le_eff"].iloc[0]) if "Le_eff" in nonlinear else np.nan
    fig.suptitle(
        rf"{100.0 * float(row['h2_volume_fraction']):.0f}% $\mathrm{{H}}_2$ nonlinear flame-front locations, "
        rf"$Le_{{\rm eff}}={le_eff:.3f}$",
        y=0.995,
        fontsize=14,
    )
    save_figure(fig, figures_dir / f"nonlinear_front_locations_{row['case']}")


def main() -> None:
    configure_reference_style()
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
        ax.plot(history["t_s"] * 1000.0, history["R_f_m"] * 1000.0, label=rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$")
    ax.set_xlabel(r"$t$ (ms)")
    ax.set_ylabel(r"$R_f$ (mm)")
    add_legend(ax, ncol=2)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Rf_t")

    plot_sb_kappa(summary, results_dir, figures_dir)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for _, row in summary.iterrows():
        history = pd.read_csv(results_dir / row["case"] / "flame_history.csv")
        ax.plot(history["R_f_m"] * 1000.0, history["S_b_m_s"], label=rf"{100 * row['h2_volume_fraction']:.0f}% $\mathrm{{H}}_2$")
    ax.set_xlabel(r"$R_f$ (mm)")
    ax.set_ylabel(r"$S_b$ (m s$^{-1}$)")
    add_legend(ax, ncol=2)
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Sb_Rf")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(summary["h2_volume_fraction"] * 100.0, summary["Ma"], "o-")
    ax.set_xlabel(r"$\mathrm{H}_2$ volume fraction (%)")
    ax.set_ylabel(r"$Ma=L_b/\delta_f$")
    ax.grid(True, alpha=0.3)
    save_figure(fig, figures_dir / "Ma_h2")

    plot_cellular_summary(summary, figures_dir)
    plot_dispersion_relations(summary, results_dir, figures_dir)
    plot_dispersion_grid_rt_modes(summary, results_dir, figures_dir)
    plot_nonlinear_summary(summary, figures_dir)
    plot_nonlinear_grid_rt_modes(summary, results_dir, figures_dir)
    plot_nonlinear_front_locations(summary, results_dir, figures_dir, config)

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
