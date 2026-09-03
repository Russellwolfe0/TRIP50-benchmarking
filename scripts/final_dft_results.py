#!/usr/bin/env python3
"""Creates a compact DFT benchmark table with reaction quantities and errors.

Run ``pivot_dft_results.py`` first. This script uses its alias-corrected
species-energy table, then writes only the useful benchmark quantities.
"""

from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIVOT_RESULTS = PROJECT_ROOT / "data/intermediate/dft_species_pivot.csv"
REFERENCE_DATA = PROJECT_ROOT / "data/intermediate/reference_data.csv"
OUTPUT_CSV = PROJECT_ROOT / "data/final/dft_benchmark_results.csv"


def number_or_none(value: str | None) -> float | None:
    return float(value) if value and value.strip() else None


def formatted(value: float | None) -> str:
    return "" if value is None else f"{value:.8f}"


def main() -> None:
    references: dict[str, dict[str, float]] = {}
    with REFERENCE_DATA.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            references[row["rxn_number"].strip()] = {
                key: float(row[key])
                for key in ("r1", "r2", "ts", "p1", "p2", "thermo", "forward_kinetic", "reverse_kinetic")
            }

    results: list[dict[str, str]] = []
    with PIVOT_RESULTS.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            reaction = row["rxn_number"].strip()
            reference = references[reaction]
            energies = {
                state: number_or_none(row.get(f"{state}_energy_kcal_mol"))
                for state in ("R1", "R2", "TS", "P1", "P2")
            }
            cpu_times = {
                state: number_or_none(row.get(f"{state}_cpu_time_s"))
                for state in ("R1", "R2", "TS", "P1", "P2")
            }

            reactants = [state for state in ("R1", "R2") if reference[state.lower()] != 0.0]
            products = [state for state in ("P1", "P2") if reference[state.lower()] != 0.0]
            required = {"TS", *reactants, *products}

            available_times = [cpu_times[state] for state in required if cpu_times[state] is not None]
            average_cpu_time = (
                sum(available_times) / len(available_times) if available_times else None
            )

            if any(energies[state] is None for state in required):
                thermo = forward = reverse = None
            else:
                reactant_energy = sum(energies[state] for state in reactants)
                product_energy = sum(energies[state] for state in products)
                thermo = product_energy - reactant_energy
                forward = energies["TS"] - reactant_energy
                reverse = energies["TS"] - product_energy

            thermo_error = thermo - reference["thermo"] if thermo is not None else None
            forward_error = forward - reference["forward_kinetic"] if forward is not None else None
            reverse_error = reverse - reference["reverse_kinetic"] if reverse is not None else None
            kinetic_error = (
                (abs(forward_error) + abs(reverse_error)) / 2
                if forward_error is not None and reverse_error is not None
                else None
            )

            results.append(
                {
                    "rxn_number": reaction,
                    "functional": row["functional"].strip(),
                    "average_cpu_time_seconds": formatted(average_cpu_time),
                    "thermo_kcal_mol": formatted(thermo),
                    "forward_kinetic_kcal_mol": formatted(forward),
                    "reverse_kinetic_kcal_mol": formatted(reverse),
                    "thermo_error_kcal_mol": formatted(thermo_error),
                    "forward_kinetic_error_kcal_mol": formatted(forward_error),
                    "reverse_kinetic_error_kcal_mol": formatted(reverse_error),
                    "kinetic_error_kcal_mol": formatted(kinetic_error),
                }
            )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
