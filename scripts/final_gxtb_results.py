#!/usr/bin/env python3
"""Create the compact gxTB benchmark table."""

from pathlib import Path

import final_mlip_results as route


PROJECT_ROOT = Path(__file__).resolve().parents[1]
route.PIVOT_RESULTS = PROJECT_ROOT / "data/intermediate/gxtb_species_pivot.csv"
route.OUTPUT_CSV = PROJECT_ROOT / "data/final/gxtb_benchmark_results.csv"
route.TIME_INPUT_LABEL = "calculation_time"
route.TIME_OUTPUT_FIELD = "average_calculation_time_seconds"


if __name__ == "__main__":
    route.main()
