#!/usr/bin/env python3
"""Plot the TRIP50 accuracy-versus-runtime Pareto front.

This is the ``final_data_table.csv`` counterpart of the archived Pareto maker:
both runtime and error are minimized, points are colored by model category, and
the nondominated points are connected in increasing-runtime order.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "data/final/final_data_table.csv"
DEFAULT_IMAGE = PROJECT_ROOT / "images/pareto_front.png"
DEFAULT_FRONT_CSV = PROJECT_ROOT / "data/final/pareto_front.csv"
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
LABEL_MODELS = {
    "aimnet2-nse": (8, 10),
    "orbmol_v1_direct": (8, -18),
    "esen-md-direct-all-omol": (8, 10),
    "B2GP-PLYP-D4": (8, -18),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create an accuracy-versus-runtime Pareto plot from the final model table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_IMAGE)
    parser.add_argument("--front-csv", type=Path, default=DEFAULT_FRONT_CSV)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    if args.dpi < 1:
        parser.error("--dpi must be >= 1")
    return args


def read_data(path: Path) -> list[dict[str, object]]:
    required = {
        "model", "model_category", "average_run_time_seconds",
        "combined_mae_kcal_mol",
    }
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        columns = set(reader.fieldnames or [])
        if not required.issubset(columns):
            raise ValueError(f"{path} is missing columns: {sorted(required - columns)}")
        rows: list[dict[str, object]] = []
        for line_number, source in enumerate(reader, start=2):
            runtime = float(source["average_run_time_seconds"])
            error = float(source["combined_mae_kcal_mol"])
            if not math.isfinite(runtime) or runtime <= 0:
                raise ValueError(f"{path}:{line_number} runtime must be finite and positive")
            if not math.isfinite(error) or error < 0:
                raise ValueError(f"{path}:{line_number} combined MAE must be finite and nonnegative")
            rows.append({
                **source,
                "runtime": runtime,
                "log_time": math.log10(runtime),
                "combined_mae": error,
            })
    if not rows:
        raise ValueError(f"No model rows found in {path}")
    return rows


def pareto_front(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Return points not strictly dominated in runtime/error space."""
    front: list[dict[str, object]] = []
    for point in rows:
        dominated = any(
            other is not point
            and float(other["runtime"]) <= float(point["runtime"])
            and float(other["combined_mae"]) <= float(point["combined_mae"])
            and (
                float(other["runtime"]) < float(point["runtime"])
                or float(other["combined_mae"]) < float(point["combined_mae"])
            )
            for other in rows
        )
        if not dominated:
            front.append(point)
    return sorted(front, key=lambda row: float(row["runtime"]))


def write_front(path: Path, front: list[dict[str, object]]) -> None:
    fields = [
        "pareto_rank", "model", "model_category", "average_run_time_seconds",
        "combined_mae_kcal_mol",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for rank, row in enumerate(front, start=1):
            writer.writerow({
                "pareto_rank": rank,
                "model": row["model"],
                "model_category": row["model_category"],
                "average_run_time_seconds": row["average_run_time_seconds"],
                "combined_mae_kcal_mol": row["combined_mae_kcal_mol"],
            })


def plot(rows: list[dict[str, object]], front: list[dict[str, object]], output: Path, dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(8, 6))
    present_categories = {str(row["model_category"]) for row in rows}
    ordered_categories = [category for category in CATEGORY_ORDER if category in present_categories]
    ordered_categories.extend(sorted(present_categories - set(ordered_categories)))
    for category in ordered_categories:
        group = [row for row in rows if row["model_category"] == category]
        axis.scatter(
            [row["log_time"] for row in group],
            [row["combined_mae"] for row in group],
            s=60,
            alpha=0.75,
            label=CATEGORY_LABELS.get(category, category),
            color=CATEGORY_COLORS.get(category, "gray"),
            edgecolors="white",
            linewidths=0.5,
            zorder=3,
        )

    axis.plot(
        [row["log_time"] for row in front],
        [row["combined_mae"] for row in front],
        color="red",
        linestyle="--",
        linewidth=1.5,
        label="Pareto front",
        zorder=2,
    )
    for row in rows:
        model = str(row["model"])
        if model not in LABEL_MODELS:
            continue
        axis.annotate(
            model,
            (float(row["log_time"]), float(row["combined_mae"])),
            xytext=LABEL_MODELS[model],
            textcoords="offset points",
            fontsize=8,
            arrowprops={"arrowstyle": "->", "linewidth": 0.7, "color": "#555555"},
        )

    axis.set_xlabel("log10(Average Run Time (s))")
    axis.set_ylabel("Combined MAE (kcal/mol)")
    axis.set_title("TRIP50 Accuracy vs Computational Cost")
    axis.grid(True, alpha=0.3)
    axis.legend()
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = read_data(args.input)
    front = pareto_front(rows)
    write_front(args.front_csv, front)
    plot(rows, front, args.output, args.dpi)
    print(f"Wrote {len(front)} Pareto-optimal models to {args.front_csv}")
    print(f"Wrote Pareto plot to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
