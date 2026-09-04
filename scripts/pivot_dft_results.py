#!/usr/bin/env python3
"""Run the shared reaction pivot for DFT CPU-time results."""

from pathlib import Path

import pivot_mlip_results as route


ROOT = Path(__file__).resolve().parents[1]
route.PARSED_RESULTS = ROOT / "data/intermediate/dft_spe_results.csv"
route.PIVOT_OUTPUT = ROOT / "data/intermediate/dft_species_pivot.csv"
route.REACTION_OUTPUT = ROOT / "data/final/dft_reaction_results.csv"
route.TIME_INPUT_FIELD = "cpu_time_seconds"
route.TIME_OUTPUT_LABEL = "cpu_time"


if __name__ == "__main__":
    route.main()
