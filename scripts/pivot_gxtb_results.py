#!/usr/bin/env python3
"""Pivot parsed gxTB results and calculate CC-reference errors."""

from pathlib import Path

import pivot_mlip_results as route


PROJECT_ROOT = Path(__file__).resolve().parents[1]
route.PARSED_RESULTS = PROJECT_ROOT / "data/intermediate/gxtb_spe_results.csv"
route.PIVOT_OUTPUT = PROJECT_ROOT / "data/intermediate/gxtb_species_pivot.csv"
route.REACTION_OUTPUT = PROJECT_ROOT / "data/final/gxtb_reaction_results.csv"
route.TIME_INPUT_FIELD = "calculation_time_seconds"
route.TIME_OUTPUT_LABEL = "calculation_time"


if __name__ == "__main__":
    route.main()
