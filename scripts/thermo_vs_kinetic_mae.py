#!/usr/bin/env python3
"""Plot model thermochemical MAE against kinetic MAE with a y=x reference."""

from __future__ import annotations

import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = PROJECT_ROOT / "data/final/final_data_table.csv"
OUTPUT_PNG = PROJECT_ROOT / "images/thermo_vs_kinetic_mae.png"
CATEGORY_ORDER = ("MLIP", "Semiempirical tight-binding", "DFT")
CATEGORY_LABELS = {
    "MLIP": "MLIP",
    "Semiempirical tight-binding": "Semi-empirical",
    "DFT": "DFT",
}
CATEGORY_COLORS = {
    "MLIP": "#2979b8",
    "Semiempirical tight-binding": "#3a923a",
    "DFT": "#7b4ab0",
}
LABEL_OFFSETS = {
    "B2GP-PLYP-D4": (8, -16),
    "wB97M-V": (8, -25),
    "aimnet2-nse": (8, 8),
    "esen-md-direct-all-omol": (8, 12),
    "orbmol_v1_direct": (8, 8),
    "gxtb": (8, -16),
}


def read_data() -> list[dict[str, object]]:
    with INPUT_CSV.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"model", "model_category", "mae_thermo_kcal_mol", "mae_kinetic_kcal_mol"}
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{INPUT_CSV} is missing columns: {sorted(required - fields)}")
        rows: list[dict[str, object]] = []
        for line_number, row in enumerate(reader, start=2):
            thermo = float(row["mae_thermo_kcal_mol"])
            kinetic = float(row["mae_kinetic_kcal_mol"])
            if not all(math.isfinite(value) and value >= 0 for value in (thermo, kinetic)):
                raise ValueError(f"{INPUT_CSV}:{line_number} contains an invalid MAE")
            rows.append({**row, "thermo": thermo, "kinetic": kinetic})
    if not rows:
        raise ValueError(f"No models found in {INPUT_CSV}")
    return rows


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = read_data()
    maximum = max(max(float(row["thermo"]), float(row["kinetic"])) for row in rows)
    upper_limit = math.ceil((maximum + 0.25) * 2) / 2

    figure, axis = plt.subplots(figsize=(8, 8))
    for category in CATEGORY_ORDER:
        group = [row for row in rows if row["model_category"] == category]
        axis.scatter(
            [row["thermo"] for row in group],
            [row["kinetic"] for row in group],
            s=72,
            alpha=0.8,
            label=CATEGORY_LABELS[category],
            color=CATEGORY_COLORS[category],
            edgecolors="white",
            linewidths=0.6,
            zorder=3,
        )

    axis.plot(
        [0, upper_limit], [0, upper_limit], color="#555555", linestyle="--",
        linewidth=1.5, label="y = x", zorder=1,
    )
    for row in rows:
        model = str(row["model"])
        if model not in LABEL_OFFSETS:
            continue
        axis.annotate(
            model,
            (float(row["thermo"]), float(row["kinetic"])),
            xytext=LABEL_OFFSETS[model],
            textcoords="offset points",
            fontsize=8,
            arrowprops={"arrowstyle": "->", "linewidth": 0.7, "color": "#555555"},
        )

    axis.set_xlim(0, upper_limit)
    axis.set_ylim(0, upper_limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Thermodynamic MAE (kcal/mol)")
    axis.set_ylabel("Kinetic MAE (kcal/mol)")
    axis.set_title("TRIP50 Thermodynamic vs Kinetic Accuracy")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Wrote MAE comparison chart to {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
