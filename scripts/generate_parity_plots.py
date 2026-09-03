#!/usr/bin/env python3
"""Generate thermodynamic and kinetic parity plots for every benchmark model."""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUTS = (
    PROJECT_ROOT / "data/final/dft_reaction_results.csv",
    PROJECT_ROOT / "data/final/mlip_reaction_results.csv",
    PROJECT_ROOT / "data/final/gxtb_reaction_results.csv",
)
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "images/parity_plots"
os.environ.setdefault("MPLBACKEND", "Agg")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate thermodynamic and kinetic parity plots per model.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("inputs", nargs="*", type=Path, help="Reaction-result CSV files")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def filename_slug(model: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", model.lower().replace("ω", "w")).strip("-")


def load_models(paths: list[Path]) -> dict[str, list[dict[str, str]]]:
    required = {
        "functional", "thermo_kcal_mol", "forward_kinetic_kcal_mol",
        "reverse_kinetic_kcal_mol", "thermo_reference_kcal_mol",
        "forward_kinetic_reference_kcal_mol", "reverse_kinetic_reference_kcal_mol",
    }
    models: dict[str, list[dict[str, str]]] = defaultdict(list)
    for path in paths:
        with path.open(newline="", encoding="utf-8-sig") as source:
            reader = csv.DictReader(source)
            fields = set(reader.fieldnames or [])
            if not required.issubset(fields):
                raise ValueError(f"{path} is missing columns: {sorted(required - fields)}")
            for row in reader:
                models[row["functional"].strip()].append(row)
    return dict(models)


def draw_plot(
    model: str,
    plot_type: str,
    reference: list[float],
    predicted: list[float],
    directions: list[str] | None,
    output_dir: Path,
) -> Path:
    import matplotlib.pyplot as plt

    lower = min(reference + predicted)
    upper = max(reference + predicted)
    padding = max((upper - lower) * 0.05, 0.5)
    limits = (lower - padding, upper + padding)
    figure, axis = plt.subplots(figsize=(6.5, 6.5), dpi=220)

    if directions is None:
        axis.scatter(
            reference, predicted, s=34, color="#1e4d2b", alpha=0.8,
            edgecolors="white", linewidths=0.4,
        )
    else:
        for direction, color in (("Forward", "#1e4d2b"), ("Reverse", "#2e86de")):
            indices = [index for index, value in enumerate(directions) if value == direction]
            axis.scatter(
                [reference[index] for index in indices],
                [predicted[index] for index in indices],
                s=34, color=color, alpha=0.8, edgecolors="white",
                linewidths=0.4, label=direction,
            )

    axis.plot(limits, limits, color="#b13a32", linestyle="--", linewidth=1.8, label="y = x")
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Reference (kcal/mol)", fontsize=12)
    axis.set_ylabel("Predicted (kcal/mol)", fontsize=12)
    axis.set_title(f"{model} {plot_type.title()} Parity Plot", fontsize=15, weight="bold")
    axis.legend(frameon=False, loc="upper left")
    axis.grid(True, color="#d9d7ca", linewidth=0.7, alpha=0.8)
    figure.tight_layout()
    output = output_dir / f"{filename_slug(model)}_{plot_type}_parity.png"
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    return output


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    paths = args.inputs or list(DEFAULT_INPUTS)
    models = load_models(paths)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expected = {
        f"{filename_slug(model)}_{plot_type}_parity.png"
        for model in models for plot_type in ("thermo", "kinetic")
    }
    for image in args.output_dir.glob("*_parity.png"):
        if image.name not in expected:
            image.unlink()

    for model in sorted(models, key=str.casefold):
        rows = models[model]
        thermo_reference = [float(row["thermo_reference_kcal_mol"]) for row in rows]
        thermo_predicted = [float(row["thermo_kcal_mol"]) for row in rows]
        draw_plot(model, "thermo", thermo_reference, thermo_predicted, None, args.output_dir)

        kinetic_reference: list[float] = []
        kinetic_predicted: list[float] = []
        directions: list[str] = []
        for row in rows:
            kinetic_reference.extend((
                float(row["forward_kinetic_reference_kcal_mol"]),
                float(row["reverse_kinetic_reference_kcal_mol"]),
            ))
            kinetic_predicted.extend((
                float(row["forward_kinetic_kcal_mol"]),
                float(row["reverse_kinetic_kcal_mol"]),
            ))
            directions.extend(("Forward", "Reverse"))
        draw_plot(
            model, "kinetic", kinetic_reference, kinetic_predicted,
            directions, args.output_dir,
        )
        print(f"Generated parity plots for {model}")

    print(f"Generated {len(models) * 2} parity plots in {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
