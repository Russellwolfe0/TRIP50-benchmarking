#!/usr/bin/env python3
"""Pivot parsed DFT single-point results and calculate CC-reference errors.

Inputs
------
``data/intermediate/dft_spe_results.csv``
    Long-form parsed output with energies in Hartree and CPU times in seconds.
``configs/aliases.csv``
    ``alias,canonical`` structure mapping. Aliases are authoritative: an alias
    gets the canonical structure's values even if it has its own calculation.
``data/intermediate/reference_data.csv``
    Coupled-cluster reference energies and reaction quantities.

Outputs
-------
``data/intermediate/dft_species_pivot.csv``
    One row per reaction, functional, and program with species energies in
    kcal/mol and per-species CPU times.
``data/final/dft_reaction_results.csv``
    Derived thermodynamic and forward/reverse kinetic quantities and errors.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


HARTREE_TO_KCAL_MOL = 627.5095
STATES = ("R1", "R2", "TS", "P1", "P2")

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PARSED_RESULTS = PROJECT_ROOT / "data/intermediate/dft_spe_results.csv"
ALIASES = PROJECT_ROOT / "configs/aliases.csv"
REFERENCE_DATA = PROJECT_ROOT / "data/intermediate/reference_data.csv"
PIVOT_OUTPUT = PROJECT_ROOT / "data/intermediate/dft_species_pivot.csv"
REACTION_OUTPUT = PROJECT_ROOT / "data/final/dft_reaction_results.csv"


def read_aliases(alias_file: Path) -> dict[str, str]:
    """Read alias -> canonical mappings, removing whitespace from old CSV rows."""
    aliases: dict[str, str] = {}
    with alias_file.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            alias = row["alias"].strip()
            canonical = row["canonical"].strip()
            if alias and canonical:
                aliases[alias] = canonical
    return aliases


def final_canonical(species: str, aliases: dict[str, str]) -> str:
    """Follow an alias chain to its canonical endpoint and reject cycles."""
    visited: set[str] = set()
    while species in aliases:
        if species in visited:
            raise ValueError(f"Alias cycle detected at {species}")
        visited.add(species)
        species = aliases[species]
    return species


def split_species(species: str) -> tuple[str, str]:
    """Split ``25-P2`` into ``('25', 'P2')``."""
    reaction, state = species.split("-", maxsplit=1)
    if state not in STATES or not reaction.isdigit():
        raise ValueError(f"Expected species in the form <reaction>-R1/R2/TS/P1/P2, got {species!r}")
    return reaction, state


def number_or_none(value: str | None) -> float | None:
    if value is None or not value.strip():
        return None
    return float(value)


def formatted(value: float | None) -> str:
    """Write missing values as blank cells and numeric values consistently."""
    return "" if value is None else f"{value:.8f}"


def state_energy(row: dict[str, str], state: str) -> float | None:
    return number_or_none(row.get(f"{state}_energy_kcal_mol"))


def reaction_quantity(row: dict[str, str], reference: dict[str, float], quantity: str) -> float | None:
    """Calculate a reaction quantity only when every physically present state exists.

    A zero energy in the reference CSV represents a state absent from that
    reaction (for example R2 in a unimolecular reaction), not a real zero-energy
    calculation. Missing states in the model data therefore remain missing.
    """
    required_reactants = [state for state in ("R1", "R2") if reference[state.lower()] != 0.0]
    required_products = [state for state in ("P1", "P2") if reference[state.lower()] != 0.0]
    required_states = {"TS", *required_reactants, *required_products}
    energies = {state: state_energy(row, state) for state in required_states}
    if any(value is None for value in energies.values()):
        return None

    reactant_energy = sum(energies[state] for state in required_reactants)
    product_energy = sum(energies[state] for state in required_products)
    transition_state_energy = energies["TS"]
    if quantity == "thermo":
        return product_energy - reactant_energy
    if quantity == "forward_kinetic":
        return transition_state_energy - reactant_energy
    if quantity == "reverse_kinetic":
        return transition_state_energy - product_energy
    raise ValueError(f"Unknown reaction quantity: {quantity}")


def main() -> None:
    for required_file in (PARSED_RESULTS, ALIASES, REFERENCE_DATA):
        if not required_file.is_file():
            raise FileNotFoundError(f"Required input file not found: {required_file}")

    aliases = read_aliases(ALIASES)

    # Store only non-alias calculations. Alias-file values are deliberately
    # ignored, even if present, because the canonical calculation takes priority.
    records: dict[tuple[str, str, str], dict[str, str]] = {}
    with PARSED_RESULTS.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            species = row["species"].strip()
            if row.get("completed", "").strip().lower() != "true" or species in aliases:
                continue
            split_species(species)
            key = (row["functional"].strip(), row["program"].strip(), species)
            if key in records:
                raise ValueError(f"Duplicate completed calculation for {key}: {records[key]['output_file']} and {row['output_file']}")
            records[key] = row

    # Build rows for every function/program pair with parsed data and every
    # reaction represented in the reference dataset.
    reference_rows: dict[str, dict[str, float]] = {}
    with REFERENCE_DATA.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            reaction = row["rxn_number"].strip()
            reference_rows[reaction] = {
                key: float(row[key])
                for key in ("r1", "r2", "ts", "p1", "p2", "thermo", "forward_kinetic", "reverse_kinetic")
            }

    function_programs = sorted({(functional, program) for functional, program, _ in records})
    pivot_rows: list[dict[str, str]] = []
    for reaction in sorted(reference_rows, key=int):
        for functional, program in function_programs:
            pivot: dict[str, str] = {
                "rxn_number": reaction,
                "functional": functional,
                "program": program,
            }
            for state in STATES:
                alias_species = f"{reaction}-{state}"
                canonical_species = final_canonical(alias_species, aliases)
                record = records.get((functional, program, canonical_species))
                energy_hartree = number_or_none(record.get("energy_hartree")) if record else None
                cpu_time = number_or_none(record.get("cpu_time_seconds")) if record else None
                pivot[f"{state}_energy_kcal_mol"] = formatted(
                    energy_hartree * HARTREE_TO_KCAL_MOL if energy_hartree is not None else None
                )
                pivot[f"{state}_cpu_time_s"] = formatted(cpu_time)
            pivot_rows.append(pivot)

    PIVOT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    pivot_fields = ["rxn_number", "functional", "program"]
    pivot_fields += [f"{state}_energy_kcal_mol" for state in STATES]
    pivot_fields += [f"{state}_cpu_time_s" for state in STATES]
    with PIVOT_OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=pivot_fields)
        writer.writeheader()
        writer.writerows(pivot_rows)

    reaction_rows: list[dict[str, str]] = []
    for pivot in pivot_rows:
        reference = reference_rows[pivot["rxn_number"]]
        thermo = reaction_quantity(pivot, reference, "thermo")
        forward = reaction_quantity(pivot, reference, "forward_kinetic")
        reverse = reaction_quantity(pivot, reference, "reverse_kinetic")
        thermo_error = thermo - reference["thermo"] if thermo is not None else None
        forward_error = forward - reference["forward_kinetic"] if forward is not None else None
        reverse_error = reverse - reference["reverse_kinetic"] if reverse is not None else None
        kinetic_error = (
            (abs(forward_error) + abs(reverse_error)) / 2
            if forward_error is not None and reverse_error is not None
            else None
        )
        reaction_rows.append(
            {
                "rxn_number": pivot["rxn_number"],
                "functional": pivot["functional"],
                "program": pivot["program"],
                "thermo_kcal_mol": formatted(thermo),
                "forward_kinetic_kcal_mol": formatted(forward),
                "reverse_kinetic_kcal_mol": formatted(reverse),
                "thermo_reference_kcal_mol": formatted(reference["thermo"]),
                "forward_kinetic_reference_kcal_mol": formatted(reference["forward_kinetic"]),
                "reverse_kinetic_reference_kcal_mol": formatted(reference["reverse_kinetic"]),
                "thermo_error_kcal_mol": formatted(thermo_error),
                "forward_kinetic_error_kcal_mol": formatted(forward_error),
                "reverse_kinetic_error_kcal_mol": formatted(reverse_error),
                "kinetic_error_kcal_mol": formatted(kinetic_error),
            }
        )

    REACTION_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    reaction_fields = list(reaction_rows[0]) if reaction_rows else []
    with REACTION_OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=reaction_fields)
        writer.writeheader()
        writer.writerows(reaction_rows)

    print(f"Wrote {len(pivot_rows)} rows to {PIVOT_OUTPUT}")
    print(f"Wrote {len(reaction_rows)} rows to {REACTION_OUTPUT}")


if __name__ == "__main__":
    main()
