#!/usr/bin/env python3
"""Run the shared benchmark aggregation for DFT CPU-time results."""

from pathlib import Path

import final_mlip_results as route


ROOT = Path(__file__).resolve().parents[1]
route.PIVOT_RESULTS = ROOT / "data/intermediate/dft_species_pivot.csv"
route.OUTPUT_CSV = ROOT / "data/final/dft_benchmark_results.csv"
route.TIME_INPUT_LABEL = "cpu_time"
route.TIME_OUTPUT_FIELD = "average_cpu_time_seconds"


if __name__ == "__main__":
    route.main()
