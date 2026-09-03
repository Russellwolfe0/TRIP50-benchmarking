#!/usr/bin/env python3
"""Generate the additional plots used by the archived-style TRIP50 website."""

from __future__ import annotations

import csv
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parents[1]
CATEGORIES = ("C-C", "C-O", "C-S", "HAT", "Si-X", "C-Hal", "N-X")
SAMPLE_REACTIONS = (1, 13, 20, 22, 31, 38, 42)
SAMPLE_MODELS = (
    "B2GP-PLYP-D4", "esen-md-direct-all-omol", "uma-m-1p1", "aimnet2-nse", "gxtb"
)
COLORS = {"DFT": "#1e4d2b", "MLIP": "#2e86de", "Semiempirical": "#b45f06"}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as source:
        return list(csv.DictReader(source))


def generate_scaled_pareto() -> None:
    import matplotlib.pyplot as plt

    data = rows(ROOT / "data/final/final_data_table.csv")
    errors = [float(row["combined_mae_kcal_mol"]) for row in data]
    low, high = min(errors), max(errors)
    scaled = [math.sqrt((value - low) / (high - low)) for value in errors]
    figure, axis = plt.subplots(figsize=(8.5, 6.3), dpi=220)
    for category, color in COLORS.items():
        selected = [(row, value) for row, value in zip(data, scaled) if row["model_category"] == category]
        axis.scatter(
            [float(row["average_run_time_seconds"]) for row, _ in selected],
            [value for _, value in selected], s=55, color=color, alpha=0.82, label=category,
            edgecolors="white", linewidths=0.5,
        )
    axis.set_xscale("log")
    axis.set_xlabel("Average calculation time (s; log scale)")
    axis.set_ylabel("Scaled combined MAE")
    axis.set_title("TRIP50 Scaled Accuracy vs Computational Cost", weight="bold")
    axis.grid(True, color="#d9d7ca", linewidth=0.7, alpha=0.8)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(ROOT / "images/paretov2_front.png", bbox_inches="tight")
    plt.close(figure)


def generate_type_radars() -> None:
    import matplotlib.pyplot as plt

    data = rows(ROOT / "data/final/reaction_type_summary.csv")
    by_category = {row["reaction_type"]: row for row in data}
    output = ROOT / "images/radar_plots"
    output.mkdir(parents=True, exist_ok=True)
    angles = [2 * math.pi * index / len(CATEGORIES) for index in range(len(CATEGORIES))]
    closed = angles + angles[:1]
    keys = {
        "DFT": "dft_mae_{kind}_kcal_mol",
        "MLIP": "mlip_mae_{kind}_kcal_mol",
        "Semiempirical": "semiempirical_mae_{kind}_kcal_mol",
    }
    for kind, title in (("thermo", "Thermodynamic"), ("kinetic", "Kinetic")):
        series = {
            model_type: [float(by_category[category][key.format(kind=kind)]) for category in CATEGORIES]
            for model_type, key in keys.items()
        }
        radial_max = max(5, math.ceil(max(map(max, series.values())) / 5) * 5)
        figure, axis = plt.subplots(figsize=(7.2, 6.7), dpi=220, subplot_kw={"polar": True})
        for model_type, values in series.items():
            axis.plot(closed, values + values[:1], color=COLORS[model_type], linewidth=2.4, label=model_type)
            axis.fill(closed, values + values[:1], color=COLORS[model_type], alpha=0.12)
        axis.set_xticks(angles)
        axis.set_xticklabels(CATEGORIES, fontsize=11)
        axis.set_ylim(0, radial_max)
        axis.set_title(f"{title} Category MAE by Method Type", fontsize=14, weight="bold", pad=22)
        axis.legend(loc="upper left", bbox_to_anchor=(1.05, 1.08), frameon=False, title="Type")
        axis.grid(color="#cfcfc6", linewidth=0.7)
        figure.tight_layout()
        figure.savefig(output / f"dft_vs_se_{kind}_categories.png", bbox_inches="tight")
        plt.close(figure)


def generate_reaction_coordinates() -> None:
    import matplotlib.pyplot as plt

    all_rows: dict[tuple[str, int], dict[str, str]] = {}
    for filename in ("dft_reaction_results.csv", "mlip_reaction_results.csv", "gxtb_reaction_results.csv"):
        for row in rows(ROOT / "data/final" / filename):
            all_rows[(row["functional"], int(row["rxn_number"]))] = row

    output_root = ROOT / "images/reaction_coordinates"
    colors = plt.get_cmap("tab10").colors
    for reaction in SAMPLE_REACTIONS:
        plotted: list[tuple[str, list[float]]] = []
        for model in SAMPLE_MODELS:
            row = all_rows[(model, reaction)]
            plotted.append((model, [0, float(row["forward_kinetic_kcal_mol"]), float(row["thermo_kcal_mol"])]))
        reference_row = all_rows[(SAMPLE_MODELS[0], reaction)]
        plotted.append(("Reference", [0, float(reference_row["forward_kinetic_reference_kcal_mol"]), float(reference_row["thermo_reference_kcal_mol"])]))
        output = output_root / f"rxn{reaction}"
        output.mkdir(parents=True, exist_ok=True)
        values = [value for _, energies in plotted for value in energies]
        figure, axis = plt.subplots(figsize=(8, 5.5), dpi=220)
        for index, (model, energies) in enumerate(plotted):
            axis.plot((0, 1, 2), energies, color=colors[index], linewidth=2.3, marker="o", markersize=5, label=model)
        axis.set_xticks((0, 1, 2), ("Reactants", "TS", "Products"))
        axis.set_ylabel("Relative energy (kcal mol$^{-1}$)")
        axis.set_title(f"Reaction {reaction}: Method Comparison")
        axis.set_xlim(-0.3, 2.3)
        axis.set_ylim(min(values) - 2, max(values) + 2)
        axis.grid(axis="y", color="#d9d7ca", linewidth=0.7)
        axis.legend(title="Method", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
        figure.tight_layout()
        figure.savefig(output / "method_comparison.png", bbox_inches="tight")
        plt.close(figure)


def main() -> int:
    generate_scaled_pareto()
    generate_type_radars()
    generate_reaction_coordinates()
    print("Generated scaled Pareto, method-type radar, and reaction-coordinate plots")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
