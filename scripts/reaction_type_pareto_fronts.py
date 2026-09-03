#!/usr/bin/env python3
"""Create thermodynamic and kinetic Pareto fronts by reaction type."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUT_CSV = PROJECT_ROOT / "data/final/model_results_by_reaction_type.csv"
OUTPUT_DIR = PROJECT_ROOT / "images/Paretofrontsbycat"
FRONT_CSV = PROJECT_ROOT / "data/final/reaction_type_pareto_fronts.csv"
ABSOLUTE_OUTPUT_DIR = PROJECT_ROOT / "images/Paretofrontsbycat_absolute"
ABSOLUTE_FRONT_CSV = PROJECT_ROOT / "data/final/reaction_type_pareto_fronts_absolute.csv"
CATEGORIES = ("C-C", "C-O", "C-S", "HAT", "Si-X", "C-Hal", "N-X")
MEASURES = {
    "thermo": ("mae_thermo_kcal_mol", "Thermodynamic MAE"),
    "kinetic": ("mae_kinetic_kcal_mol", "Kinetic MAE"),
}
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


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create reaction-type thermodynamic and kinetic Pareto plots.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--absolute", action="store_true",
        help="Plot absolute MAE in kcal/mol and write to separate absolute-MAE outputs",
    )
    return parser.parse_args(argv)


def read_data() -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    with INPUT_CSV.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {
            "reaction_type", "model", "model_category", "average_run_time_seconds",
            "mae_thermo_kcal_mol", "mae_kinetic_kcal_mol",
        }
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{INPUT_CSV} is missing columns: {sorted(required - fields)}")
        for line_number, row in enumerate(reader, start=2):
            runtime = float(row["average_run_time_seconds"])
            if not math.isfinite(runtime) or runtime <= 0:
                raise ValueError(f"{INPUT_CSV}:{line_number} has invalid runtime")
            grouped[row["reaction_type"]].append({
                **row,
                "runtime": runtime,
                "log_time": math.log10(runtime),
                "thermo": float(row["mae_thermo_kcal_mol"]),
                "kinetic": float(row["mae_kinetic_kcal_mol"]),
            })
    for category in CATEGORIES:
        if len(grouped[category]) != 24:
            raise ValueError(f"Expected 24 models for {category}, found {len(grouped[category])}")
    return grouped


def pareto_front(rows: list[dict[str, object]], measure: str) -> list[dict[str, object]]:
    front: list[dict[str, object]] = []
    for point in rows:
        dominated = any(
            other is not point
            and float(other["runtime"]) <= float(point["runtime"])
            and float(other[measure]) <= float(point[measure])
            and (
                float(other["runtime"]) < float(point["runtime"])
                or float(other[measure]) < float(point[measure])
            )
            for other in rows
        )
        if not dominated:
            front.append(point)
    return sorted(front, key=lambda row: float(row["runtime"]))


def plot_front(
    rows: list[dict[str, object]], front: list[dict[str, object]],
    measure: str, label: str, category: str, scaled: bool, output_dir: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    raw_values = [float(row[measure]) for row in rows]
    minimum = min(raw_values)
    value_range = max(raw_values) - minimum
    value_key = "plot_mae"
    for row in rows:
        row[value_key] = (
            100 * (float(row[measure]) - minimum) / value_range if value_range else 0.0
        ) if scaled else float(row[measure])

    figure, axis = plt.subplots(figsize=(8, 6))
    for model_category in ("MLIP", "Semiempirical tight-binding", "DFT"):
        group = [row for row in rows if row["model_category"] == model_category]
        axis.scatter(
            [row["log_time"] for row in group],
            [row[value_key] for row in group],
            s=60,
            alpha=0.75,
            label=CATEGORY_LABELS[model_category],
            color=CATEGORY_COLORS[model_category],
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )
    axis.plot(
        [row["log_time"] for row in front],
        [row[value_key] for row in front],
        color="red", linestyle="--", linewidth=1.5, zorder=2, label="Pareto front",
    )
    label_y: dict[str, float] = {}
    previous: float | None = None
    label_gap = 5.5 if scaled else max(value_range * 0.055, 0.15)
    for row in sorted(front, key=lambda item: float(item[value_key])):
        point_y = float(row[value_key])
        placed_y = point_y if previous is None else max(point_y, previous + label_gap)
        label_y[str(row["model"])] = placed_y
        previous = placed_y
    for row in front:
        point_x = float(row["log_time"])
        place_left = point_x > 2.0
        axis.annotate(
            str(row["model"]),
            (point_x, float(row[value_key])),
            xytext=(point_x - 0.12 if place_left else point_x + 0.12, label_y[str(row["model"])]),
            textcoords="data",
            fontsize=7,
            ha="right" if place_left else "left",
            va="center",
            arrowprops={"arrowstyle": "->", "linewidth": 0.7, "color": "#555555"},
        )
    axis.set_xlabel(f"log10(Average time for {category} reactions (s))")
    axis.set_ylabel(
        f"Scaled {label} (% of model range)" if scaled else f"{label} (kcal/mol)"
    )
    axis.set_title(
        f"TRIP50 Scaled {label}: {category}" if scaled else f"TRIP50 {label}: {category}"
    )
    axis.legend()
    axis.grid(True, alpha=0.3)
    figure.tight_layout()
    slug = re.sub(r"[^a-z0-9]+", "_", category.lower()).strip("_")
    output = output_dir / f"pareto_{measure}_{slug}.png"
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    grouped = read_data()
    output_dir = ABSOLUTE_OUTPUT_DIR if args.absolute else OUTPUT_DIR
    front_csv = ABSOLUTE_FRONT_CSV if args.absolute else FRONT_CSV
    output_dir.mkdir(parents=True, exist_ok=True)
    front_rows: list[dict[str, object]] = []
    for measure, (_, label) in MEASURES.items():
        for category in CATEGORIES:
            rows = grouped[category]
            front = pareto_front(rows, measure)
            plot_front(rows, front, measure, label, category, not args.absolute, output_dir)
            for rank, row in enumerate(front, start=1):
                front_rows.append({
                    "reaction_type": category,
                    "measure": measure,
                    "pareto_rank": rank,
                    "model": row["model"],
                    "model_category": row["model_category"],
                    "average_run_time_seconds": f"{float(row['runtime']):.8f}",
                    "mae_kcal_mol": f"{float(row[measure]):.8f}",
                })
    fields = [
        "reaction_type", "measure", "pareto_rank", "model", "model_category",
        "average_run_time_seconds", "mae_kcal_mol",
    ]
    front_csv.parent.mkdir(parents=True, exist_ok=True)
    with front_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(front_rows)
    print(f"Wrote {len(CATEGORIES) * len(MEASURES)} plots to {output_dir}")
    print(f"Wrote {len(front_rows)} Pareto-front rows to {front_csv}")


if __name__ == "__main__":
    main()
