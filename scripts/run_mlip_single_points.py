#!/usr/bin/env python3
"""Run consistent ASE single-point calculations with every TRIP50 MLIP."""
"""Example CLI: python scripts/run_mlip_single_points.py \
  --backend mace \
  --model extra_large \
  --output calculations/mlip_results/mace-omol-extra-large.csv"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
import time
from pathlib import Path
from typing import Any

from mlip_calculator_factories import BACKENDS, build_provider

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = PROJECT_ROOT / "configs/structure_metadata.csv"
EV_TO_HARTREE = 1.0 / 27.211386245988
XYZ_SUFFIXES = {".xyz"}
FIELDS = [
    "species", "model", "backend", "energy_eV", "energy_hartree",
    "max_force_eV_per_A", "n_atoms", "charge", "spin_multiplicity", "device",
    "repeats", "mean_time_seconds", "median_time_seconds", "min_time_seconds",
    "stdev_time_seconds", "status", "input_file", "error",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Universal ASE energy/force driver for TRIP50 MLIPs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--backend", required=True, choices=sorted(BACKENDS))
    parser.add_argument("--model", required=True, help="Checkpoint/model identifier")
    parser.add_argument("--input", type=Path, default=PROJECT_ROOT / "structures/reference")
    parser.add_argument(
        "--metadata", type=Path, default=DEFAULT_METADATA,
        help="CSV containing species, charge, and multiplicity",
    )
    parser.add_argument("--output", type=Path, help="CSV path (derived from model if omitted)")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--charge", type=int, help="Override CSV charge for every structure")
    parser.add_argument(
        "--multiplicity", type=int,
        help="Override CSV multiplicity for every structure",
    )
    parser.add_argument("--no-recursive", action="store_true")
    parser.add_argument("--forces", action="store_true", help="Evaluate forces as well as energy")
    parser.add_argument("--warmup", type=int, default=1)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument(
        "--calculator-kwargs", default="{}", metavar="JSON",
        help='Backend options, e.g. \'{"family":"polar","default_dtype":"float64"}\'',
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.warmup < 0 or args.repeats < 1:
        parser.error("--warmup must be >= 0 and --repeats must be >= 1")
    if args.multiplicity is not None and args.multiplicity < 1:
        parser.error("--multiplicity must be >= 1")
    try:
        args.calculator_kwargs = json.loads(args.calculator_kwargs)
    except json.JSONDecodeError as error:
        parser.error(f"invalid --calculator-kwargs JSON: {error}")
    if not isinstance(args.calculator_kwargs, dict):
        parser.error("--calculator-kwargs must decode to a JSON object")
    return args


def files_to_run(path: Path, recursive: bool = True) -> list[Path]:
    if path.is_file():
        return [path] if path.suffix.lower() in XYZ_SUFFIXES else []
    if not path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    candidates = path.rglob("*") if recursive else path.glob("*")
    return sorted(p for p in candidates if p.is_file() and p.suffix.lower() in XYZ_SUFFIXES)


def load_structure_metadata(path: Path) -> dict[str, tuple[int, int]]:
    """Read and validate species -> (charge, multiplicity) metadata."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Structure metadata not found: {path}. Run scripts/chargemult_parser.py first."
        )

    metadata: dict[str, tuple[int, int]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"species", "charge", "multiplicity"}
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{path} is missing columns: {sorted(required - fields)}")

        for row_number, row in enumerate(reader, start=2):
            species = row["species"].strip()
            if not species:
                raise ValueError(f"{path}:{row_number} has an empty species")
            if species in metadata:
                raise ValueError(f"{path}:{row_number} duplicates species {species!r}")
            try:
                charge = int(row["charge"])
                multiplicity = int(row["multiplicity"])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{row_number} has a non-integer charge or multiplicity"
                ) from error
            if multiplicity < 1:
                raise ValueError(f"{path}:{row_number} has multiplicity {multiplicity}; expected >= 1")
            metadata[species] = (charge, multiplicity)

    if not metadata:
        raise ValueError(f"No structure metadata rows found in {path}")
    return metadata


def synchronize(device: str) -> None:
    if not device.startswith("cuda"):
        return
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.synchronize()
    except ImportError:
        pass


def load_atoms(
    path: Path,
    metadata: dict[str, tuple[int, int]],
    charge_override: int | None = None,
    multiplicity_override: int | None = None,
):
    try:
        from ase.io import read
    except ImportError as error:
        raise RuntimeError("ASE is not installed in the active model environment") from error
    try:
        metadata_charge, metadata_multiplicity = metadata[path.stem]
    except KeyError as error:
        raise KeyError(
            f"No metadata for species {path.stem!r} in the selected metadata CSV"
        ) from error
    atoms = read(path)
    atoms.info["charge"] = metadata_charge if charge_override is None else charge_override
    atoms.info["spin"] = (
        metadata_multiplicity if multiplicity_override is None else multiplicity_override
    )
    return atoms


def evaluate(atoms: Any, calculator: Any, device: str, forces: bool, repeats: int) -> dict[str, Any]:
    timings: list[float] = []
    energy = 0.0
    max_force: float | str = ""
    atoms.calc = calculator
    for _ in range(repeats):
        if hasattr(calculator, "reset"):
            calculator.reset()
        synchronize(device)
        started = time.perf_counter()
        if forces:
            force_array = atoms.get_forces()
            energy = float(atoms.get_potential_energy())
            max_force = float((force_array**2).sum(axis=1).max() ** 0.5)
        else:
            energy = float(atoms.get_potential_energy())
        synchronize(device)
        timings.append(time.perf_counter() - started)
    return {
        "energy_eV": f"{energy:.12f}",
        "energy_hartree": f"{energy * EV_TO_HARTREE:.12f}",
        "max_force_eV_per_A": "" if max_force == "" else f"{max_force:.12f}",
        "mean_time_seconds": f"{statistics.mean(timings):.8f}",
        "median_time_seconds": f"{statistics.median(timings):.8f}",
        "min_time_seconds": f"{min(timings):.8f}",
        "stdev_time_seconds": f"{statistics.stdev(timings) if len(timings) > 1 else 0.0:.8f}",
    }


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(PROJECT_ROOT.resolve()))
    except ValueError:
        return str(path)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    files = files_to_run(args.input, not args.no_recursive)
    if not files:
        raise FileNotFoundError(f"No .xyz or .extxyz files found at {args.input}")
    metadata = load_structure_metadata(args.metadata)
    unlisted = [path for path in files if path.stem not in metadata]
    for path in unlisted:
        print(f"Skipping {path}: no row in {args.metadata}", file=sys.stderr)
    files = [path for path in files if path.stem in metadata]
    if not files:
        raise ValueError(f"None of the selected structures have rows in {args.metadata}")
    output = args.output or PROJECT_ROOT / "data/intermediate/mlip_results" / f"{args.model}.csv"
    print(f"{args.backend}/{args.model}: {len(files)} structure(s) on {args.device}")
    if args.dry_run:
        for path in files:
            charge, multiplicity = metadata[path.stem]
            print(f"{path} | charge={charge} multiplicity={multiplicity}")
        return 0

    provider = build_provider(args.backend, args.model, args.device, **args.calculator_kwargs)
    for path in files[:args.warmup]:
        atoms = load_atoms(path, metadata, args.charge, args.multiplicity)
        evaluate(atoms, provider.for_atoms(atoms), args.device, args.forces, 1)

    rows: list[dict[str, Any]] = []
    for index, path in enumerate(files, 1):
        atoms = None
        row = {field: "" for field in FIELDS}
        row.update({
            "species": path.stem, "model": args.model, "backend": args.backend,
            "device": args.device, "repeats": args.repeats, "status": "failed",
            "input_file": display_path(path),
        })
        try:
            atoms = load_atoms(path, metadata, args.charge, args.multiplicity)
            row.update({
                "n_atoms": len(atoms), "charge": atoms.info["charge"],
                "spin_multiplicity": atoms.info["spin"],
            })
            row.update(evaluate(
                atoms, provider.for_atoms(atoms), args.device, args.forces, args.repeats
            ))
            row["status"] = "completed"
            print(f"[{index}/{len(files)}] {path.name}: {row['energy_eV']} eV")
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
            print(f"[{index}/{len(files)}] FAILED {path.name}: {row['error']}", file=sys.stderr)
        rows.append(row)

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    completed = sum(row["status"] == "completed" for row in rows)
    print(f"Wrote {len(rows)} rows ({completed} completed) to {output}")
    return 0 if completed == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
