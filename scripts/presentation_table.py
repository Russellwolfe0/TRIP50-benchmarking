#!/usr/bin/env python3
"""Render the complete TRIP50 model summary as a presentation-ready PNG."""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = PROJECT_ROOT / "data/final/final_data_table.csv"
OUTPUT_PNG = PROJECT_ROOT / "images/final_data_table.png"


def display_runtime(seconds: float) -> str:
    if seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with INPUT_CSV.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 24:
        raise ValueError(f"Expected 24 models in {INPUT_CSV}, found {len(rows)}")
    rows.sort(key=lambda row: float(row["combined_mae_kcal_mol"]))

    category_names = {
        "DFT": "DFT",
        "MLIP": "MLIP",
        "Semiempirical tight-binding": "Semi-empirical",
    }
    headers = [
        "Model", "Category", "Avg. runtime", "Thermo MAE",
        "Kinetic MAE", "Combined MAE",
    ]
    table_rows: list[list[str]] = []
    for row in rows:
        table_rows.append([
            row["model"],
            category_names.get(row["model_category"], row["model_category"]),
            display_runtime(float(row["average_run_time_seconds"])),
            f"{float(row['mae_thermo_kcal_mol']):.2f}",
            f"{float(row['mae_kinetic_kcal_mol']):.2f}",
            f"{float(row['combined_mae_kcal_mol']):.2f}",
        ])

    figure, axis = plt.subplots(figsize=(16, 7.2))
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    axis.axis("off")

    table = axis.table(
        cellText=table_rows,
        colLabels=headers,
        cellLoc="left",
        colLoc="left",
        colWidths=[0.30, 0.14, 0.14, 0.13, 0.14, 0.15],
        bbox=[0.01, 0.01, 0.98, 0.98],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)

    category_tints = {
        "DFT": "#f0eafa",
        "MLIP": "#e8f2fa",
        "Semi-empirical": "#eaf5e9",
    }
    category_text = {
        "DFT": "#69409a",
        "MLIP": "#216a9a",
        "Semi-empirical": "#347a3b",
    }
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.PAD = 0.025
        if row_index == 0:
            cell.set_facecolor("#263442")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
            cell.get_text().set_fontsize(10.5)
            continue
        category = table_rows[row_index - 1][1]
        base = category_tints.get(category, "#ffffff")
        cell.set_facecolor(base if row_index % 2 else "#ffffff")
        cell.get_text().set_color("#202b36")
        if column_index == 0:
            cell.get_text().set_fontweight("semibold")
        elif column_index == 1:
            cell.get_text().set_color(category_text.get(category, "#536171"))
            cell.get_text().set_fontweight("bold")
        elif column_index >= 2:
            cell.get_text().set_ha("right")
        if column_index == 5:
            cell.get_text().set_fontweight("bold")

    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(
        OUTPUT_PNG, dpi=240, facecolor=figure.get_facecolor(),
        bbox_inches="tight", pad_inches=0.03,
    )
    plt.close(figure)
    print(f"Wrote presentation table to {OUTPUT_PNG}")


if __name__ == "__main__":
    main()
