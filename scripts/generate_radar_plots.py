#!/usr/bin/env python3
"""Generate thermodynamic and kinetic reaction-category radar plots per model.

The current reaction-category table is stored in long form: one row per model
and reaction type. This script pivots those rows in memory and writes one PNG
per model for use by the TRIP50 website.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "data/final/model_results_by_reaction_type.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "images/radar_plots"
CATEGORIES = ("C-C", "C-O", "C-S", "HAT", "Si-X", "C-Hal", "N-X")

# Render PNGs on local machines, clusters, and CI without a display server.
os.environ.setdefault("MPLBACKEND", "Agg")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate one reaction-category MAE radar plot per model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def filename_slug(model: str) -> str:
    """Create stable, browser-safe filenames while retaining readable names."""
    return re.sub(r"[^a-z0-9]+", "-", model.lower().replace("ω", "w")).strip("-")


def load_models(path: Path) -> dict[str, dict[str, dict[str, float]]]:
    """Load and validate the long-form reaction-category result table."""
    required = {
        "reaction_type",
        "model",
        "mae_thermo_kcal_mol",
        "mae_kinetic_kcal_mol",
    }
    models: dict[str, dict[str, dict[str, float]]] = {}
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{path} is missing columns: {sorted(required - fields)}")
        for line_number, row in enumerate(reader, start=2):
            model = row["model"].strip()
            category = row["reaction_type"].strip()
            if category not in CATEGORIES:
                raise ValueError(f"{path}:{line_number} has unknown reaction type {category!r}")
            values = models.setdefault(model, {"thermo": {}, "kinetic": {}})
            if category in values["thermo"]:
                raise ValueError(f"{path}:{line_number} duplicates {model!r}/{category!r}")
            values["thermo"][category] = float(row["mae_thermo_kcal_mol"])
            values["kinetic"][category] = float(row["mae_kinetic_kcal_mol"])

    for model, properties in models.items():
        for property_name, values in properties.items():
            missing = set(CATEGORIES) - set(values)
            if missing:
                raise ValueError(
                    f"{model!r} is missing {property_name} values for {sorted(missing)}"
                )
    return models


def main(argv: list[str] | None = None) -> int:
    import matplotlib.pyplot as plt

    args = parse_args(argv)
    models = load_models(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    expected_files = {f"radar_{filename_slug(model)}.png" for model in models}
    for image_file in args.output_dir.glob("radar_*.png"):
        if image_file.name not in expected_files:
            image_file.unlink()

    angles = [2 * math.pi * index / len(CATEGORIES) for index in range(len(CATEGORIES))]
    closed_angles = angles + angles[:1]

    for model in sorted(models, key=str.casefold):
        thermo = [models[model]["thermo"][category] for category in CATEGORIES]
        kinetic = [models[model]["kinetic"][category] for category in CATEGORIES]
        radial_max = max(5, math.ceil(max(thermo + kinetic) / 5) * 5)

        figure, axis = plt.subplots(
            figsize=(7.2, 6.7), dpi=220, subplot_kw={"polar": True}
        )
        axis.plot(
            closed_angles, thermo + thermo[:1], color="#1e4d2b",
            linewidth=2.2, label="Thermodynamic MAE",
        )
        axis.fill(closed_angles, thermo + thermo[:1], color="#1e4d2b", alpha=0.14)
        axis.plot(
            closed_angles, kinetic + kinetic[:1], color="#2e86de",
            linewidth=2.2, label="Kinetic MAE",
        )
        axis.fill(closed_angles, kinetic + kinetic[:1], color="#2e86de", alpha=0.14)
        axis.set_xticks(angles)
        axis.set_xticklabels(CATEGORIES, fontsize=11)
        axis.set_ylim(0, radial_max)
        axis.set_yticks(list(range(5, radial_max + 1, 5)))
        axis.set_yticklabels(
            [str(value) for value in range(5, radial_max + 1, 5)], fontsize=8
        )
        axis.set_title(
            f"{model} Reaction-Category MAE", fontsize=14, weight="bold", pad=22
        )
        axis.legend(loc="upper left", bbox_to_anchor=(1.08, 1.1), frameon=False)
        axis.grid(color="#cfcfc6", linewidth=0.7)
        figure.tight_layout()
        output = args.output_dir / f"radar_{filename_slug(model)}.png"
        figure.savefig(output, bbox_inches="tight")
        plt.close(figure)
        print(f"Wrote {output}")

    print(f"Generated {len(models)} radar plots in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
