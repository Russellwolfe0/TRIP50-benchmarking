#!/usr/bin/env python3
"""Normalize gxTB driver results for the TRIP50 benchmark pipeline."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "calculations/gxtb_results/gxtb.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "data/intermediate/gxtb_spe_results.csv"
FIELDS = [
    "species", "functional", "program", "energy_hartree",
    "calculation_time_seconds", "completed", "output_file",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert gxTB driver CSVs into normalized species results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs", nargs="*", type=Path,
        help="gxTB result CSV file(s) or directories searched recursively",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def csv_files(inputs: Iterable[Path]) -> list[Path]:
    files: set[Path] = set()
    for path in inputs:
        if path.is_file():
            files.add(path.resolve())
        elif path.is_dir():
            files.update(candidate.resolve() for candidate in path.rglob("*.csv"))
        else:
            raise FileNotFoundError(f"gxTB input does not exist: {path}")
    return sorted(files)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = args.inputs or [DEFAULT_INPUT]
    files = [path for path in csv_files(inputs) if path != args.output.resolve()]
    if not files:
        raise FileNotFoundError(f"No gxTB CSV files found under: {', '.join(map(str, inputs))}")

    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()
    for path in files:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            required = {"species", "energy_hartree", "mean_time_seconds", "status"}
            fields = set(reader.fieldnames or [])
            if not required.issubset(fields):
                raise ValueError(f"{path} is missing columns: {sorted(required - fields)}")
            count = 0
            for line_number, row in enumerate(reader, start=2):
                species = row["species"].strip()
                functional = (row.get("model") or "gxtb").strip()
                if not species:
                    raise ValueError(f"{path}:{line_number} has no species")
                key = (functional, species)
                if key in seen:
                    raise ValueError(f"Duplicate gxTB result for {functional}/{species}")
                seen.add(key)
                rows.append({
                    "species": species,
                    "functional": functional,
                    "program": "gxTB",
                    "energy_hartree": row["energy_hartree"].strip(),
                    "calculation_time_seconds": row["mean_time_seconds"].strip(),
                    "completed": str(row["status"].strip().lower() == "completed"),
                    "output_file": row.get("output_file", "").strip(),
                })
                count += 1
        print(f"Parsed {count} rows from {path}")

    rows.sort(key=lambda row: (row["functional"], row["species"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
