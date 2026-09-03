#!/usr/bin/env python3
"""Combine standardized MLIP result CSVs for downstream benchmarking.

Input files are produced by ``run_mlip_single_points.py``. The output columns
The runner synchronizes the GPU immediately before and after each calculation,
so ``mean_time_seconds`` is carried through as GPU calculation time rather than
being mislabeled as CPU or general wall time.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "calculations/mlip_results"
DEFAULT_OUTPUT = PROJECT_ROOT / "data/intermediate/mlip_spe_results.csv"
FIELDS = [
    "species",
    "functional",
    "program",
    "energy_hartree",
    "wall_time_seconds",
    "gpu_time_seconds",
    "cpu_time_seconds",
    "cpu_time_source",
    "nprocs",
    "completed",
    "output_file",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert MLIP result CSVs to the dft_spe_results.csv schema.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        help="MLIP CSV file(s) or directories; directories are searched recursively",
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
            raise FileNotFoundError(f"MLIP input does not exist: {path}")
    return sorted(files)


def first_value(row: dict[str, str], *names: str) -> str:
    for name in names:
        value = row.get(name, "").strip()
        if value:
            return value
    return ""


def completed_value(row: dict[str, str], energy: str) -> bool:
    status = row.get("status", "").strip().lower()
    if status:
        return status in {"completed", "complete", "success", "ok", "true"}
    completed = row.get("completed", "").strip().lower()
    if completed:
        return completed in {"true", "1", "yes"}
    return bool(energy)


def relative_or_original(path_text: str) -> str:
    if not path_text:
        return ""
    path = Path(path_text)
    if not path.is_absolute():
        return path.as_posix()
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_file(path: Path) -> list[dict[str, object]]:
    results: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            return results
        required = {"species", "energy_hartree"}
        if not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} is not a standardized MLIP result CSV; missing "
                f"{sorted(required - set(reader.fieldnames))}"
            )
        for row_number, row in enumerate(reader, start=2):
            species = row.get("species", "").strip()
            model = first_value(row, "model", "functional") or path.stem
            energy = row.get("energy_hartree", "").strip()
            if not species:
                raise ValueError(f"{path}:{row_number} has no species value")
            results.append(
                {
                    "species": species,
                    "functional": model,
                    "program": "ASE",
                    "energy_hartree": energy,
                    "wall_time_seconds": first_value(row, "wall_time_seconds"),
                    "gpu_time_seconds": first_value(
                        row, "gpu_time_seconds", "mean_time_seconds", "time_seconds", "time_s"
                    ),
                    "cpu_time_seconds": "",
                    "cpu_time_source": "",
                    "nprocs": "",
                    "completed": completed_value(row, energy),
                    "output_file": relative_or_original(row.get("input_file", "").strip()),
                }
            )
    return results


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    inputs = args.inputs or [DEFAULT_INPUT]
    files = csv_files(inputs)
    output_resolved = args.output.resolve()
    files = [path for path in files if path != output_resolved]
    if not files:
        raise FileNotFoundError(f"No MLIP CSV files found under: {', '.join(map(str, inputs))}")

    rows: list[dict[str, object]] = []
    for path in files:
        parsed = parse_file(path)
        rows.extend(parsed)
        print(f"Parsed {len(parsed)} rows from {path}")
    rows.sort(key=lambda row: (str(row["functional"]), str(row["species"])))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
