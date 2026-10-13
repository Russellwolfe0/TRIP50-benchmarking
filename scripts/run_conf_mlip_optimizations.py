#!/usr/bin/env python3
"""Geometry-optimize every AQME/RDKit conformer of every TRIP50 species with every MLIP model.

Mirrors the backend abstraction used by TRIP50's single-point MLIP driver
(``mlip_calculator_factories.py``, copied unmodified from
``TRIP50/scripts/``): each backend package (aimnet, fairchem, mace, orb) is
imported lazily. Since no single environment is expected to have all four
backend packages installed, pass the model(s) available in whatever
environment you currently have activated with ``--models`` (or edit
SELECTED_MODELS below), then repeat with different models in each other
model environment.

Inputs (all in geom/, the folder above scripts/):
    conformer_xyz/<species>.xyz   multi-frame XYZ, one frame per conformer,
                                  made from reference_conformers/<species>_rdkit.sdf with
                                  ``obabel <sdf> -O conformer_xyz/<species>.xyz`` (no -m).
                                  Only files named exactly <species>.xyz are read, so
                                  stray per-conformer files (e.g. 3-R11.xyz from -m) are ignored.
    structure_metadata.csv        charge and multiplicity for every species. These are always
                                  used instead of the SDF <Mult> tag, which is wrong (singlet)
                                  for the 13 excited-state triplets.
    mlip_calculator_factories.py  calculator adapters.

Transition states have no conformers (a SMILES can't define a TS), so they're skipped.

For each selected model, every conformer frame gets a BFGS optimization, writing
    conformer_optimizations/<model>/<species>_conf<i>_optimized.xyz
    conformer_optimizations/<model>/<species>_conf<i>.traj / .log
and one CSV per model, conformer_optimization_results/<model>.csv, written
incrementally so a crash mid-run only costs the conformers not yet written.
``conformer_index`` is the 0-based frame index in conformer_xyz/<species>.xyz
(not RDKit's conformer ID shown in the XYZ comment line). ``energy_eV`` is the
final energy, used later to pick the lowest-energy conformer per species per model.

Timing excludes kernel warmup: the first conformer run per model gets one
throwaway energy evaluation before the timer starts, so one-time
kernel/graph compilation isn't counted against it.

Reruns resume automatically: any (species, conformer_index) pair already
marked "completed" in a model's results CSV is skipped.

Data lives in calculations/conformer_study/ (override with TRIP50_STUDY_DIR).

For now this only generates results -- no aggregation or analysis.

Run from anywhere, e.g. on ladanyi with a thread cap:
    OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 nohup nice -n 10 \\
        python3 run_conf_mlip_optimizations.py --models esen-md-direct-all-omol uma-m-1p1 \\
        > conf_fairchem.log 2>&1 &

Quick test on a couple of species first:
    python3 run_conf_mlip_optimizations.py --models esen-md-direct-all-omol --species 1-R1 3-R1
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import os
from pathlib import Path
from typing import Any

REPO_DIR = Path(__file__).resolve().parents[1]  # scripts/ lives at the TRIP50 repo root
# Conformer-study data (raw runs and derived tables). Override with TRIP50_STUDY_DIR to point at
# a different copy, e.g. the large raw optimizations kept outside the repo.
STUDY_DIR = Path(os.environ.get("TRIP50_STUDY_DIR", REPO_DIR / "calculations" / "conformer_study"))

from mlip_calculator_factories import build_provider  # sits next to this script in scripts/

CONFORMER_DIR = STUDY_DIR / "conformer_xyz"
METADATA_CSV = REPO_DIR / "configs" / "structure_metadata.csv"
OUTPUT_ROOT = STUDY_DIR / "conformer_optimizations"
RESULTS_DIR = STUDY_DIR / "conformer_optimization_results"

DEVICE = "cuda"
FMAX = 0.05
MAX_STEPS = 500

# model -> backend for all 16 TRIP50 MLIPs, matching
# TRIP50/calculations/mlip_results/*.csv, grouped by backend.
MODELS = {
    "aimnet2-nse": "aimnet",

    "allscaip-md-conserving-all-omol": "fairchem",
    "allscaip-md-direct-all-omol": "fairchem",
    "esen-md-direct-all-omol": "fairchem",
    "esen-sm-conserving-all-omol": "fairchem",
    "esen-sm-direct-all-omol": "fairchem",
    "uma-m-1p1": "fairchem",
    "uma-s-1p1": "fairchem",
    "uma-s-1p2": "fairchem",

    "MACE-OMOL": "mace",
    "polar-1-l": "mace",
    "polar-1-m": "mace",
    "polar-1-s": "mace",

    "orbmol_v1_conservative": "orb",
    "orbmol_v1_direct": "orb",
    "orbmol_v2": "orb",
}

# Default when --models isn't given. Running with all 16 selected will fail on
# any model whose backend isn't installed in the current environment.
SELECTED_MODELS = list(MODELS)

FIELDS = [
    "species", "conformer_index", "model", "backend", "n_atoms",
    "charge", "spin_multiplicity", "device", "n_steps", "converged",
    "fmax_final", "energy_eV", "elapsed_time_seconds", "status", "input_file",
    "output_file", "error",
]


def load_metadata(path: Path) -> dict[str, tuple[int, int]]:
    metadata = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            metadata[row["species"]] = (int(row["charge"]), int(row["multiplicity"]))
    return metadata


def conformers_to_optimize(
    metadata: dict[str, tuple[int, int]],
    only_species: set[str] | None = None,
) -> list[tuple[str, str, Path]]:
    """Return (species, conformer_index, xyz_path) for every conformer frame.

    conformer_index is a string so it matches the value read back from the
    results CSV when resuming.
    """
    from ase.io import read

    jobs = []
    for species in sorted(metadata):
        if only_species and species not in only_species:
            continue
        path = CONFORMER_DIR / f"{species}.xyz"
        if not path.exists():  # TS species have no conformers
            continue
        n_frames = len(read(path, index=":"))
        jobs.extend((species, str(index), path) for index in range(n_frames))
    return jobs


def load_completed(results_csv: Path) -> dict[tuple[str, str], dict[str, Any]]:
    """Return {(species, conformer_index): row} for previously completed runs."""
    if not results_csv.exists():
        return {}
    with results_csv.open(newline="", encoding="utf-8") as stream:
        return {
            (row["species"], row["conformer_index"]): row
            for row in csv.DictReader(stream)
            if row.get("status") == "completed"
        }


def synchronize(device: str) -> None:
    if not device.startswith("cuda"):
        return
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.synchronize()
    except ImportError:
        pass


def run_model(
    model: str,
    backend: str,
    provider: Any,
    metadata: dict[str, tuple[int, int]],
    jobs: list[tuple[str, str, Path]],
    device: str,
    writer: csv.DictWriter,
    stream: Any,
) -> int:
    """Run every conformer for one model, writing (and flushing) each row as it completes."""

    def emit(row: dict[str, Any]) -> None:
        writer.writerow(row)
        stream.flush()

    from ase.io import read, write
    from ase.optimize import BFGS

    output_dir = OUTPUT_ROOT / model
    output_dir.mkdir(parents=True, exist_ok=True)

    completed = 0
    warmed_up = False
    for index, (species, conf_index, xyz_file) in enumerate(jobs, 1):
        charge, multiplicity = metadata[species]
        tag = f"{species}_conf{conf_index}"
        row = {field: "" for field in FIELDS}
        row.update({
            "species": species, "conformer_index": conf_index, "model": model,
            "backend": backend, "charge": charge, "spin_multiplicity": multiplicity,
            "device": device, "status": "failed", "input_file": str(xyz_file),
        })

        try:
            atoms = read(xyz_file, index=int(conf_index))
            atoms.info["charge"] = charge
            atoms.info["spin"] = multiplicity
            atoms.calc = provider.for_atoms(atoms)

            if not warmed_up:
                atoms.get_potential_energy()
                synchronize(device)
                if hasattr(atoms.calc, "reset"):
                    atoms.calc.reset()
                warmed_up = True

            optimizer = BFGS(
                atoms,
                trajectory=str(output_dir / f"{tag}.traj"),
                logfile=str(output_dir / f"{tag}.log"),
            )
            synchronize(device)
            started = time.perf_counter()
            optimizer.run(fmax=FMAX, steps=MAX_STEPS)
            synchronize(device)
            elapsed = time.perf_counter() - started

            fmax_final = float((atoms.get_forces() ** 2).sum(axis=1).max() ** 0.5)
            energy = float(atoms.get_potential_energy())
            output_file = output_dir / f"{tag}_optimized.xyz"
            write(output_file, atoms)

            row.update({
                "n_atoms": len(atoms),
                "n_steps": optimizer.get_number_of_steps(),
                "converged": fmax_final <= FMAX,
                "fmax_final": f"{fmax_final:.6f}",
                "energy_eV": f"{energy:.8f}",
                "elapsed_time_seconds": f"{elapsed:.4f}",
                "status": "completed",
                "output_file": str(output_file),
            })
            completed += 1
            print(f"[{model}] [{index}/{len(jobs)}] {tag}: {elapsed:.2f}s, {row['n_steps']} steps, E = {energy:.6f} eV")
        except Exception as error:
            row["error"] = f"{type(error).__name__}: {error}"
            print(f"[{model}] [{index}/{len(jobs)}] FAILED {tag}: {row['error']}", file=sys.stderr)

        emit(row)

    return completed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--models", nargs="+", default=SELECTED_MODELS, choices=list(MODELS), metavar="MODEL",
        help="model(s) to run; must be installed in the current environment (default: all 16)",
    )
    parser.add_argument(
        "--species", nargs="+", default=None, metavar="SPECIES",
        help="only run these species, e.g. --species 1-R1 3-R1 (for testing)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    metadata = load_metadata(METADATA_CSV)
    jobs = conformers_to_optimize(metadata, set(args.species) if args.species else None)
    if not jobs:
        raise FileNotFoundError(f"No conformer structures found to optimize in {CONFORMER_DIR}")

    n_species = len({job[0] for job in jobs})
    print(f"{len(jobs)} conformer(s) across {n_species} species")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    for model in args.models:
        backend = MODELS[model]
        results_csv = RESULTS_DIR / f"{model}.csv"
        already_done = load_completed(results_csv)
        remaining_jobs = [job for job in jobs if (job[0], job[1]) not in already_done]

        skipped = len(jobs) - len(remaining_jobs)
        note = f" ({skipped} already completed, skipping)" if skipped else ""
        print(f"=== {model} ({backend}): {len(remaining_jobs)} conformer(s) to run{note} on {DEVICE} ===")
        if not remaining_jobs:
            print(f"{model}: nothing to do, all {len(jobs)} already completed")
            continue

        try:
            provider = build_provider(backend, model, DEVICE)
        except Exception as error:
            print(
                f"Failed to build provider for {model} ({backend}): "
                f"{type(error).__name__}: {error}",
                file=sys.stderr,
            )
            continue

        with results_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader()
            for row in already_done.values():
                writer.writerow(row)
            stream.flush()
            completed = run_model(model, backend, provider, metadata, remaining_jobs, DEVICE, writer, stream)
        print(f"{model}: {completed}/{len(remaining_jobs)} newly completed -> {results_csv}")


if __name__ == "__main__":
    main()
