#!/usr/bin/env python3
"""Calculate and render TRIP50 metrics separated by reaction category."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REACTION_TYPES = PROJECT_ROOT / "configs/rxntypes.csv"
REFERENCE_DATA = PROJECT_ROOT / "data/intermediate/reference_data.csv"
MODEL_OUTPUT = PROJECT_ROOT / "data/final/model_results_by_reaction_type.csv"
CATEGORY_OUTPUT = PROJECT_ROOT / "data/final/reaction_type_summary.csv"
IMAGE_DIR = PROJECT_ROOT / "images/reaction_types"
CATEGORY_IMAGE = PROJECT_ROOT / "images/reaction_type_summary.png"
STATES = ("R1", "R2", "TS", "P1", "P2")
ROUTES = (
    (
        PROJECT_ROOT / "data/intermediate/dft_species_pivot.csv",
        "DFT", "cpu_time",
    ),
    (
        PROJECT_ROOT / "data/intermediate/mlip_species_pivot.csv",
        "MLIP", "gpu_time",
    ),
    (
        PROJECT_ROOT / "data/intermediate/gxtb_species_pivot.csv",
        "Semiempirical tight-binding", "calculation_time",
    ),
)
MODEL_FIELDS = [
    "reaction_type", "model", "model_category", "average_run_time_seconds",
    "mae_thermo_kcal_mol", "mae_kinetic_kcal_mol", "combined_mae_kcal_mol",
]
CATEGORY_FIELDS = [
    "reaction_type", "mae_thermo_kcal_mol", "mae_kinetic_kcal_mol",
    "dft_mae_thermo_kcal_mol", "dft_mae_kinetic_kcal_mol",
    "mlip_mae_thermo_kcal_mol", "mlip_mae_kinetic_kcal_mol",
    "semiempirical_mae_thermo_kcal_mol", "semiempirical_mae_kinetic_kcal_mol",
    "average_run_time_seconds",
]


def number(value: str | None) -> float | None:
    return float(value) if value and value.strip() else None


def mean(values: list[float]) -> float:
    if not values:
        raise ValueError("Cannot average an empty collection")
    return sum(values) / len(values)


def formatted(value: float) -> str:
    return f"{value:.8f}"


def read_reaction_types() -> dict[str, str]:
    mapping: dict[str, str] = {}
    with REACTION_TYPES.open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        next(reader, None)
        for line_number, row in enumerate(reader, start=2):
            if len(row) < 2:
                raise ValueError(f"{REACTION_TYPES}:{line_number} has fewer than two columns")
            reaction = row[0].strip()
            reaction_type = row[1].strip()
            if not reaction or not reaction_type:
                raise ValueError(f"{REACTION_TYPES}:{line_number} has an empty value")
            if reaction in mapping:
                raise ValueError(f"{REACTION_TYPES}:{line_number} duplicates reaction {reaction}")
            mapping[reaction] = reaction_type
    return mapping


def read_references() -> dict[str, dict[str, float]]:
    references: dict[str, dict[str, float]] = {}
    with REFERENCE_DATA.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            references[row["rxn_number"].strip()] = {
                key: float(row[key])
                for key in (
                    "r1", "r2", "ts", "p1", "p2", "thermo",
                    "forward_kinetic", "reverse_kinetic",
                )
            }
    return references


def reaction_metrics(
    row: dict[str, str], reference: dict[str, float], time_label: str,
) -> tuple[float, float, float] | None:
    reactants = [state for state in ("R1", "R2") if reference[state.lower()] != 0.0]
    products = [state for state in ("P1", "P2") if reference[state.lower()] != 0.0]
    required = {"TS", *reactants, *products}
    energies = {
        state: number(row.get(f"{state}_energy_kcal_mol")) for state in required
    }
    times = {state: number(row.get(f"{state}_{time_label}_s")) for state in required}
    if any(value is None for value in energies.values()):
        return None
    available_times = [value for value in times.values() if value is not None]
    if not available_times:
        return None
    reactant_energy = sum(float(energies[state]) for state in reactants)
    product_energy = sum(float(energies[state]) for state in products)
    transition_energy = float(energies["TS"])
    thermo_error = (product_energy - reactant_energy) - reference["thermo"]
    forward_error = (transition_energy - reactant_energy) - reference["forward_kinetic"]
    reverse_error = (transition_energy - product_energy) - reference["reverse_kinetic"]
    kinetic_error = (abs(forward_error) + abs(reverse_error)) / 2
    return abs(thermo_error), kinetic_error, mean([float(value) for value in available_times])


def calculate() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    reaction_types = read_reaction_types()
    references = read_references()
    if set(references) != set(reaction_types):
        missing_types = sorted(set(references) - set(reaction_types), key=int)
        missing_references = sorted(set(reaction_types) - set(references), key=int)
        raise ValueError(
            f"Reaction mapping mismatch; missing types={missing_types}, "
            f"missing references={missing_references}"
        )

    grouped: dict[tuple[str, str, str], list[tuple[float, float, float]]] = defaultdict(list)
    for pivot_path, model_category, time_label in ROUTES:
        if not pivot_path.is_file():
            raise FileNotFoundError(f"Required pivot not found: {pivot_path}")
        with pivot_path.open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                reaction = row["rxn_number"].strip()
                metrics = reaction_metrics(row, references[reaction], time_label)
                if metrics is not None:
                    key = (reaction_types[reaction], row["functional"].strip(), model_category)
                    grouped[key].append(metrics)

    model_rows: list[dict[str, str]] = []
    for (reaction_type, model, model_category), values in grouped.items():
        thermo_mae = mean([value[0] for value in values])
        kinetic_mae = mean([value[1] for value in values])
        run_time = mean([value[2] for value in values])
        model_rows.append({
            "reaction_type": reaction_type,
            "model": model,
            "model_category": model_category,
            "average_run_time_seconds": formatted(run_time),
            "mae_thermo_kcal_mol": formatted(thermo_mae),
            "mae_kinetic_kcal_mol": formatted(kinetic_mae),
            "combined_mae_kcal_mol": formatted((thermo_mae + kinetic_mae) / 2),
        })
    model_rows.sort(key=lambda row: (row["reaction_type"], float(row["combined_mae_kcal_mol"])))

    by_category: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in model_rows:
        by_category[row["reaction_type"]].append(row)
    category_rows: list[dict[str, str]] = []
    for reaction_type, rows in sorted(by_category.items()):
        model_categories = {
            "dft": [row for row in rows if row["model_category"] == "DFT"],
            "mlip": [row for row in rows if row["model_category"] == "MLIP"],
            "semiempirical": [
                row for row in rows
                if row["model_category"] == "Semiempirical tight-binding"
            ],
        }
        category_maes: dict[str, str] = {}
        for prefix, selected in model_categories.items():
            if not selected:
                raise ValueError(f"No {prefix} models found for reaction type {reaction_type}")
            category_maes[f"{prefix}_mae_thermo_kcal_mol"] = formatted(
                mean([float(row["mae_thermo_kcal_mol"]) for row in selected])
            )
            category_maes[f"{prefix}_mae_kinetic_kcal_mol"] = formatted(
                mean([float(row["mae_kinetic_kcal_mol"]) for row in selected])
            )
        category_rows.append({
            "reaction_type": reaction_type,
            "mae_thermo_kcal_mol": formatted(mean([float(row["mae_thermo_kcal_mol"]) for row in rows])),
            "mae_kinetic_kcal_mol": formatted(mean([float(row["mae_kinetic_kcal_mol"]) for row in rows])),
            **category_maes,
            "average_run_time_seconds": formatted(mean([float(row["average_run_time_seconds"]) for row in rows])),
        })
    return model_rows, category_rows


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def display_runtime(seconds: float) -> str:
    if seconds < 1:
        return f"{seconds * 1000:.1f} ms"
    if seconds < 120:
        return f"{seconds:.1f} s"
    if seconds < 7200:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def render_table(
    rows: list[list[str]], headers: list[str], output: Path,
    column_widths: list[float], category_column: int | None = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    height = max(2.2, 0.31 * (len(rows) + 1))
    figure, axis = plt.subplots(figsize=(16, height))
    figure.patch.set_facecolor("white")
    axis.axis("off")
    table = axis.table(
        cellText=rows, colLabels=headers, cellLoc="left", colLoc="left",
        colWidths=column_widths, bbox=[0.01, 0.01, 0.98, 0.98],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.5)
    category_colors = {
        "DFT": ("#f0eafa", "#69409a"),
        "MLIP": ("#e8f2fa", "#216a9a"),
        "Semi-empirical": ("#eaf5e9", "#347a3b"),
    }
    for (row_index, column_index), cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.PAD = 0.025
        if row_index == 0:
            cell.set_facecolor("#263442")
            cell.get_text().set_color("white")
            cell.get_text().set_fontweight("bold")
            continue
        category = rows[row_index - 1][category_column] if category_column is not None else ""
        tint, text_color = category_colors.get(category, ("#edf1f4", "#202b36"))
        cell.set_facecolor(tint if row_index % 2 else "white")
        cell.get_text().set_color("#202b36")
        if column_index == 0:
            cell.get_text().set_fontweight("semibold")
        if category_column is not None and column_index == category_column:
            cell.get_text().set_color(text_color)
            cell.get_text().set_fontweight("bold")
        if column_index >= (2 if category_column is not None else 1):
            cell.get_text().set_ha("right")
        if column_index == len(headers) - 1 and category_column is not None:
            cell.get_text().set_fontweight("bold")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=240, bbox_inches="tight", pad_inches=0.03)
    plt.close(figure)


def render(model_rows: list[dict[str, str]], category_rows: list[dict[str, str]]) -> None:
    category_names = {
        "DFT": "DFT", "MLIP": "MLIP",
        "Semiempirical tight-binding": "Semi-empirical",
    }
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in model_rows:
        grouped[row["reaction_type"]].append(row)
    for reaction_type, rows in grouped.items():
        display_rows = [[
            row["model"], category_names[row["model_category"]],
            display_runtime(float(row["average_run_time_seconds"])),
            f"{float(row['mae_thermo_kcal_mol']):.2f}",
            f"{float(row['mae_kinetic_kcal_mol']):.2f}",
            f"{float(row['combined_mae_kcal_mol']):.2f}",
        ] for row in rows]
        filename = re.sub(r"[^a-z0-9]+", "_", reaction_type.lower()).strip("_") + ".png"
        render_table(
            display_rows,
            ["Model", "Category", "Avg. runtime", "Thermo MAE", "Kinetic MAE", "Combined MAE"],
            IMAGE_DIR / filename,
            [0.30, 0.14, 0.14, 0.13, 0.14, 0.15],
            category_column=1,
        )

    summary_rows = [[
        row["reaction_type"],
        f"{float(row['mae_thermo_kcal_mol']):.2f}",
        f"{float(row['dft_mae_thermo_kcal_mol']):.2f}",
        f"{float(row['mlip_mae_thermo_kcal_mol']):.2f}",
        f"{float(row['semiempirical_mae_thermo_kcal_mol']):.2f}",
        f"{float(row['mae_kinetic_kcal_mol']):.2f}",
        f"{float(row['dft_mae_kinetic_kcal_mol']):.2f}",
        f"{float(row['mlip_mae_kinetic_kcal_mol']):.2f}",
        f"{float(row['semiempirical_mae_kinetic_kcal_mol']):.2f}",
        display_runtime(float(row["average_run_time_seconds"])),
    ] for row in category_rows]
    render_table(
        summary_rows,
        [
            "Reaction type", "Thermo all", "Thermo DFT", "Thermo MLIP", "Thermo SE",
            "Kinetic all", "Kinetic DFT", "Kinetic MLIP", "Kinetic SE", "Avg. runtime",
        ],
        CATEGORY_IMAGE,
        [0.14, 0.095, 0.095, 0.095, 0.095, 0.095, 0.095, 0.095, 0.095, 0.10],
    )


def main() -> None:
    model_rows, category_rows = calculate()
    write_csv(MODEL_OUTPUT, MODEL_FIELDS, model_rows)
    write_csv(CATEGORY_OUTPUT, CATEGORY_FIELDS, category_rows)
    render(model_rows, category_rows)
    print(f"Wrote {len(model_rows)} model/category rows to {MODEL_OUTPUT}")
    print(f"Wrote {len(category_rows)} reaction-category rows to {CATEGORY_OUTPUT}")
    print(f"Wrote category tables to {IMAGE_DIR}")
    print(f"Wrote category summary table to {CATEGORY_IMAGE}")


if __name__ == "__main__":
    main()
