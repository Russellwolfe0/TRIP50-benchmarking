#!/usr/bin/env python3
"""Build authoritative TRIP50 charge/multiplicity metadata from updated XYZ files."""

from __future__ import annotations

from pathlib import Path
import csv
import re

PROJECT_ROOT = Path(__file__).resolve().parents[1]

XYZ_FILES = PROJECT_ROOT / "structures/reference"
OUTPUT_CSV = PROJECT_ROOT / "configs/structure_metadata.csv"
ALIASES_CSV = PROJECT_ROOT / "configs/aliases.csv"
REFERENCE_CSV = PROJECT_ROOT / "data/intermediate/reference_data.csv"

STATES = ("R1", "R2", "TS", "P1", "P2")


def read_aliases(path: Path) -> dict[str, str]:
    aliases = {}

    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            alias = row["alias"].strip()
            canonical = row["canonical"].strip()

            if alias and canonical:
                aliases[alias] = canonical
    
    return aliases


def final_canonical(
        species: str,
        aliases: dict[str, str],
) -> str:
    visited = set()

    while species in aliases:
        if species in visited:
            raise ValueError(f"Aliases cycle detected at {species}")
        
        visited.add(species)
        species = aliases[species]

    return species


def required_species(
        reference_csv: Path,
        aliases: dict[str, str],
) -> set[str]:
    required = set()

    with reference_csv.open(newline="", encoding = "utf-8") as stream:
        for row in csv.DictReader(stream):
            reaction = row["rxn_number"].strip()

            for state in STATES:
                if float(row[state.lower()]) != 0.0:
                    species = f"{reaction}-{state}"
                    required.add(final_canonical(species, aliases))
                                 
    return required


FULL_LABELS = {
    "singlet": 1,
    "doublet": 2,
    "triplet": 3,
}

ABBREVIATIONS = {
    "sing": 1,
    "doub": 2,
    "trip": 3,
}


def xyz_comment(path: Path) -> str:
    with path.open(encoding="utf-8", errors="replace") as stream:
        atom_count = next(stream, "").strip()
        comment = next(stream, "").strip()

    if not atom_count:
        raise ValueError(f"{path} is empty")

    try:
        int(atom_count)
    except ValueError as error:
        raise ValueError(
            f"{path}: line 1 is not an atom count"
        ) from error

    if not comment:
        raise ValueError(f"{path}: line 2 has no comment")

    return comment


def multiplicity_from_comment(
    comment: str,
    source: Path,
) -> int:
    normalized = comment.lower()
    matches = set()

    for label, multiplicity in FULL_LABELS.items():
        # Full electronic-state names are unambiguous even when attached to an
        # underscore or trailing zero (for example ``_triplet``/``triplet0``).
        if label in normalized:
            matches.add(multiplicity)

    for label, multiplicity in ABBREVIATIONS.items():
        pattern = rf"(?:^|[-_]){label}(?:[-_0-9]|$)"

        if re.search(pattern, normalized):
            matches.add(multiplicity)

    if not matches:
        raise ValueError(
            f"{source}: could not determine multiplicity "
            f"from comment {comment!r}"
        )

    if len(matches) > 1:
        raise ValueError(
            f"{source}: conflicting multiplicity labels "
            f"in comment {comment!r}: {sorted(matches)}"
        )

    return next(iter(matches))


def main() -> None:
    for required_path in (XYZ_FILES, ALIASES_CSV, REFERENCE_CSV):
        if not required_path.exists():
            raise FileNotFoundError(f"Required input not found: {required_path}")

    aliases = read_aliases(ALIASES_CSV)
    required = required_species(REFERENCE_CSV, aliases)
    if len(required) != 156:
        raise ValueError(
            f"Expected 156 canonical TRIP50 structures, found {len(required)}"
        )

    rows: list[dict[str, object]] = []
    for species in sorted(required):
        xyz_file = XYZ_FILES / f"{species}.xyz"
        if not xyz_file.is_file():
            raise FileNotFoundError(f"Missing required structure: {xyz_file}")

        comment = xyz_comment(xyz_file)
        rows.append(
            {
                "species": species,
                # Every TRIP50 calculation is neutral.
                "charge": 0,
                "multiplicity": multiplicity_from_comment(comment, xyz_file),
                "xyz_file": str(xyz_file.relative_to(PROJECT_ROOT)),
                "comment": comment,
            }
        )

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species", "charge", "multiplicity", "xyz_file", "comment"]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
