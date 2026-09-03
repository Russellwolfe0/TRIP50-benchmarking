#!/usr/bin/env python3
"""Run organized g-xTB single-point calculations for TRIP50 structures.

Each structure gets an isolated raw-output directory containing its input,
charge/spin control files, and complete g-xTB output.``.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PROJECT_ROOT / "structures/reference"
DEFAULT_METADATA = PROJECT_ROOT / "configs/structure_metadata.csv"
DEFAULT_WORK_DIR = PROJECT_ROOT / "calculations/gxtb/raw_outputs"
DEFAULT_OUTPUT = PROJECT_ROOT / "calculations/gxtb_results/gxtb.csv"
HARTREE_TO_EV = 27.211386245988
ENERGY_RE = re.compile(r"^\s*total\s+(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)\s*$", re.MULTILINE)
FIELDS = [
    "species", "model", "backend", "energy_eV", "energy_hartree",
    "n_atoms", "charge", "spin_multiplicity", "device", "mean_time_seconds",
    "status", "input_file", "output_file", "error",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run g-xTB single points for TRIP50 XYZ structures.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--timeout", type=float, help="Per-calculation timeout in seconds")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    if args.timeout is not None and args.timeout <= 0:
        parser.error("--timeout must be > 0")
    return args


def xyz_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path] if path.suffix.lower() == ".xyz" else []
    if not path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    candidates = path.rglob("*.xyz")
    return sorted(candidate for candidate in candidates if candidate.is_file())


def load_metadata(path: Path) -> dict[str, tuple[int, int]]:
    metadata: dict[str, tuple[int, int]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"species", "charge", "multiplicity"}
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            raise ValueError(f"{path} is missing columns: {sorted(required - fields)}")
        for line_number, row in enumerate(reader, start=2):
            species = row["species"].strip()
            if not species:
                raise ValueError(f"{path}:{line_number} has an empty species")
            if species in metadata:
                raise ValueError(f"{path}:{line_number} duplicates {species!r}")
            charge = int(row["charge"])
            multiplicity = int(row["multiplicity"])
            if multiplicity < 1:
                raise ValueError(f"{path}:{line_number} has multiplicity < 1")
            metadata[species] = (charge, multiplicity)
    return metadata


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def atom_count(path: Path) -> int:
    with path.open(encoding="utf-8") as stream:
        try:
            return int(stream.readline().strip())
        except ValueError as error:
            raise ValueError(f"Invalid XYZ atom count in {path}") from error


def prepare_directory(source: Path, run_dir: Path, charge: int, multiplicity: int) -> Path:
    run_dir.mkdir(parents=True, exist_ok=True)
    local_xyz = run_dir / source.name
    shutil.copy2(source, local_xyz)
    (run_dir / ".CHRG").write_text(f"{charge}\n", encoding="utf-8")
    (run_dir / ".UHF").write_text(f"{multiplicity - 1}\n", encoding="utf-8")
    return local_xyz


def run_once(
    executable: str, source: Path, run_dir: Path, charge: int, multiplicity: int,
    timeout: float | None,
) -> tuple[float, float]:
    local_xyz = prepare_directory(source, run_dir, charge, multiplicity)
    command = [executable, "-c", local_xyz.name]
    command.append("--gxtb")
    started = time.perf_counter()
    try:
        result = subprocess.run(
            command, cwd=run_dir, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=timeout, check=False,
        )
    except subprocess.TimeoutExpired as error:
        output = error.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        (run_dir / "gxtb.out").write_text(output, encoding="utf-8")
        raise RuntimeError(f"g-xTB timed out after {timeout:g} seconds") from error
    elapsed = time.perf_counter() - started
    (run_dir / "gxtb.out").write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(f"g-xTB exited with status {result.returncode}")
    matches = ENERGY_RE.findall(result.stdout)
    if not matches:
        raise RuntimeError("Could not locate final total energy in gxtb.out")
    return float(matches[-1]), elapsed


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    files = xyz_files(args.input)
    if not files:
        raise FileNotFoundError(f"No XYZ files found at {args.input}")
    metadata = load_metadata(DEFAULT_METADATA)
    unlisted = [path for path in files if path.stem not in metadata]
    for path in unlisted:
        print(f"Skipping {path}: no row in {DEFAULT_METADATA}", file=sys.stderr)
    files = [path for path in files if path.stem in metadata]
    if not files:
        raise ValueError(f"None of the selected structures have rows in {DEFAULT_METADATA}")

    executable = "gxtb"
    if not args.dry_run and shutil.which(executable) is None:
        raise FileNotFoundError("g-xTB executable not found on PATH: 'gxtb'")

    print(f"gxtb: {len(files)} structure(s)")
    if args.dry_run:
        for path in files:
            charge, multiplicity = metadata[path.stem]
            print(f"{path} | charge={charge} multiplicity={multiplicity}")
        return 0

    rows: list[dict[str, object]] = []
    for index, path in enumerate(files, start=1):
        charge, multiplicity = metadata[path.stem]
        species_dir = args.work_dir / path.stem
        error_text = ""
        energy_hartree = None
        elapsed = None
        try:
            energy_hartree, elapsed = run_once(
                executable, path, species_dir, charge, multiplicity,
                args.timeout,
            )
        except Exception as error:
            error_text = f"{type(error).__name__}: {error}"

        completed = energy_hartree is not None
        row = {field: "" for field in FIELDS}
        row.update({
            "species": path.stem, "model": "gxtb", "backend": "gxtb",
            "energy_eV": "" if energy_hartree is None else f"{energy_hartree * HARTREE_TO_EV:.12f}",
            "energy_hartree": "" if energy_hartree is None else f"{energy_hartree:.12f}",
            "n_atoms": atom_count(path), "charge": charge,
            "spin_multiplicity": multiplicity, "device": "cpu",
            "mean_time_seconds": "" if elapsed is None else f"{elapsed:.8f}",
            "status": "completed" if completed else "failed",
            "input_file": display_path(path),
            "output_file": display_path(species_dir / "gxtb.out"), "error": error_text,
        })
        rows.append(row)
        detail = f"{row['energy_hartree']} Eh" if completed else error_text
        print(f"[{index}/{len(files)}] {path.name}: {detail}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    completed_count = sum(row["status"] == "completed" for row in rows)
    print(f"Wrote {len(rows)} rows ({completed_count} completed) to {args.output}")
    return 0 if completed_count == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
