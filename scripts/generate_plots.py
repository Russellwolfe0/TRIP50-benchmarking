#!/usr/bin/env python3
"""Generate TRIP50 figures from one command."""

from __future__ import annotations

import argparse

import generate_parity_plots
import generate_radar_plots
import generate_website_plots
import pareto_front
import presentation_table
import reaction_category_tables
import reaction_type_pareto_fronts
import thermo_vs_kinetic_mae


TASKS = (
    "tables", "pareto", "reaction-pareto", "reaction-pareto-absolute",
    "accuracy", "parity", "radar", "website",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate TRIP50 tables and plots.")
    parser.add_argument("tasks", nargs="*", choices=TASKS, help="default: all")
    return parser.parse_args()


def main() -> None:
    selected = set(parse_args().tasks or TASKS)
    if "tables" in selected:
        reaction_category_tables.main()
        presentation_table.main()
    if "pareto" in selected:
        pareto_front.main([])
    if "reaction-pareto" in selected:
        reaction_type_pareto_fronts.main([])
    if "reaction-pareto-absolute" in selected:
        reaction_type_pareto_fronts.main(["--absolute"])
    if "accuracy" in selected:
        thermo_vs_kinetic_mae.main()
    if "parity" in selected:
        generate_parity_plots.main([])
    if "radar" in selected:
        generate_radar_plots.main([])
    if "website" in selected:
        generate_website_plots.main()


if __name__ == "__main__":
    main()
