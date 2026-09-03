#!/usr/bin/env python3
"""Aggregate all TRIP50 benchmark routes into one model-level data table.

For each model, runtime is averaged over reaction-level average calculation
times. Thermochemical MAE is the mean absolute thermochemical error, kinetic
MAE is the mean of the per-reaction forward/reverse kinetic MAE, and combined
MAE gives the thermochemical and kinetic MAEs equal weight.
"""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FINAL_DIR = PROJECT_ROOT / "data/final"
OUTPUT_CSV = FINAL_DIR / "final_data_table.csv"
ROUTES = (
    (FINAL_DIR / "dft_benchmark_results.csv", "DFT", "average_cpu_time_seconds"),
    (FINAL_DIR / "mlip_benchmark_results.csv", "MLIP", "average_gpu_time_seconds"),
    (
        FINAL_DIR / "gxtb_benchmark_results.csv",
        "Semiempirical tight-binding",
        "average_calculation_time_seconds",
    ),
)
FIELDS = [
    "model",
    "model_category",
    "average_run_time_seconds",
    "mae_thermo_kcal_mol",
    "mae_kinetic_kcal_mol",
    "combined_mae_kcal_mol",
]


def numeric(value: str | None) -> float | None:
    return float(value) if value and value.strip() else None


def mean(values: list[float], model: str, quantity: str) -> float:
    if not values:
        raise ValueError(f"No {quantity} values available for {model}")
    return sum(values) / len(values)


def formatted(value: float) -> str:
    return f"{value:.8f}"


def main() -> None:
    summaries: list[dict[str, str]] = []
    seen_models: set[str] = set()

    for input_csv, category, time_field in ROUTES:
        if not input_csv.is_file():
            raise FileNotFoundError(f"Required benchmark table not found: {input_csv}")
        grouped: dict[str, list[dict[str, str]]] = {}
        with input_csv.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            required = {
                "functional", time_field, "thermo_error_kcal_mol",
                "kinetic_error_kcal_mol",
            }
            columns = set(reader.fieldnames or [])
            if not required.issubset(columns):
                raise ValueError(f"{input_csv} is missing columns: {sorted(required - columns)}")
            for row in reader:
                grouped.setdefault(row["functional"].strip(), []).append(row)

        for model, rows in grouped.items():
            if model in seen_models:
                raise ValueError(f"Model {model!r} occurs in more than one benchmark category")
            seen_models.add(model)
            run_times = [value for row in rows if (value := numeric(row[time_field])) is not None]
            thermo_errors = [
                abs(value)
                for row in rows
                if (value := numeric(row["thermo_error_kcal_mol"])) is not None
            ]
            kinetic_errors = [
                abs(value)
                for row in rows
                if (value := numeric(row["kinetic_error_kcal_mol"])) is not None
            ]
            average_run_time = mean(run_times, model, "runtime")
            mae_thermo = mean(thermo_errors, model, "thermochemical error")
            mae_kinetic = mean(kinetic_errors, model, "kinetic error")
            combined_mae = (mae_thermo + mae_kinetic) / 2
            summaries.append({
                "model": model,
                "model_category": category,
                "average_run_time_seconds": formatted(average_run_time),
                "mae_thermo_kcal_mol": formatted(mae_thermo),
                "mae_kinetic_kcal_mol": formatted(mae_kinetic),
                "combined_mae_kcal_mol": formatted(combined_mae),
            })

    category_order = {"DFT": 0, "MLIP": 1, "Semiempirical tight-binding": 2}
    summaries.sort(key=lambda row: (category_order[row["model_category"]], row["model"].lower()))
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(summaries)
    print(f"Wrote {len(summaries)} model summaries to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
