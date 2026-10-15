#!/usr/bin/env python3
"""Analyze the MLIP conformer optimizations against the TRIP50 reference structures.

Reads the optimization *output files* rather than the per-model results CSVs
(some of those CSVs have fused rows from overlapping runs), so it only needs:

    (paths below are relative to STUDY_DIR = calculations/conformer_study unless noted)
    conformer_optimizations/<model>/<species>_conf<i>_optimized.xyz   final geometry + energy= tag (eV)
    conformer_optimizations/<model>/<species>_conf<i>.log             ASE BFGS log (steps, final fmax)
    optimizations/<model>/<species>_optimized.xyz                     same model, optimized from the
                                                                       reference geometry (distortion study)
    <repo>/structures/reference/<species>.xyz                         DFT reference geometry
    <repo>/configs/structure_metadata.csv                             charge / multiplicity

For every model x species x conformer it records the energy, convergence,
heavy-atom RMSD to the DFT reference (atoms matched by graph, via RDKit
GetBestRMS on perceived connectivity, so differing atom order is fine), and
whether the bonding pattern (including H) is still the same as the reference.

For every model x species it takes the lowest-energy optimized conformer and
compares it with the optimized reference from the same model:

    dE_vs_ref_kcal = E(best conformer) - E(optimized reference)

and assigns one category (first match wins):

    rearranged          bonding (incl. H) differs from the reference
    not_converged       best conformer did not reach fmax <= 0.05 eV/A
    linear_artifact     best conformer has more near-linear (>170 deg) two-coordinate
                        heavy-atom centres than the reference -- the UFF/MMFF starting
                        geometry trapped BFGS on a symmetric saddle
    matches_reference   |dE| <= 0.5 kcal/mol
    lower_than_reference  dE < -0.5 kcal/mol   (candidate new conformer, or MLIP artefact)
    higher_than_reference dE > +0.5 kcal/mol   (reference conformer not recovered)

Outputs (in conformer_analysis/):
    conformer_table.csv   one row per model x species x conformer
    best_conformers.csv   one row per model x species (lowest-energy conformer + category)
    model_summary.csv     one row per model
    species_consensus.csv one row per species: how the models agree

Run (in the aqme env, which has RDKit):
    python analyze_conformer_optimizations.py
    python analyze_conformer_optimizations.py --models uma-m-1p1 uma-s-1p1
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from collections import Counter, defaultdict
import os
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolAlign

RDLogger.DisableLog("rdApp.*")

REPO_DIR = Path(__file__).resolve().parents[1]  # scripts/ lives at the TRIP50 repo root
# Conformer-study data (raw runs and derived tables). Override with TRIP50_STUDY_DIR to point at
# a different copy, e.g. the large raw optimizations kept outside the repo.
STUDY_DIR = Path(os.environ.get("TRIP50_STUDY_DIR", REPO_DIR / "calculations" / "conformer_study"))
CONF_OPT_DIR = STUDY_DIR / "conformer_optimizations"
REF_OPT_DIR = STUDY_DIR / "optimizations"
REFERENCE_DIR = REPO_DIR / "structures" / "reference"
METADATA_CSV = REPO_DIR / "configs" / "structure_metadata.csv"
OUTPUT_DIR = STUDY_DIR / "conformer_analysis"

EV_TO_KCAL = 23.060548
FMAX = 0.05
DE_TOL = 0.5          # kcal/mol window counted as "matches reference"
LINEAR_ANGLE = 170.0  # degrees
MAX_MATCHES = 20000   # symmetry-equivalent mappings tried by GetBestRMS

COVALENT_RADII = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
                  "Si": 1.11, "S": 1.05, "Cl": 1.02, "Br": 1.20}

CONF_FILE = re.compile(r"^(?P<species>\d+-(?:R1|R2|P1|P2))_conf(?P<index>\d+)_optimized\.xyz$")


# ---------- reading ----------

def read_xyz(path: Path) -> tuple[list[str], np.ndarray, float | None]:
    """Return (elements, coordinates, energy_eV or None) from a (extended) XYZ file."""
    lines = path.read_text().splitlines()
    n_atoms = int(lines[0])
    match = re.search(r"energy=(-?[\d.eE+-]+)", lines[1])
    elements, coords = [], []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        elements.append(parts[0])
        coords.append([float(v) for v in parts[1:4]])
    return elements, np.array(coords), (float(match.group(1)) if match else None)


def read_bfgs_log(path: Path) -> tuple[int | None, float | None]:
    """Return (steps, final fmax) from the last line of an ASE BFGS log."""
    if not path.exists():
        return None, None
    last = [line for line in path.read_text().splitlines() if line.startswith("BFGS:")]
    if not last:
        return None, None
    parts = last[-1].split()
    return int(parts[1]), float(parts[-1])


def load_metadata() -> dict[str, tuple[int, int]]:
    with METADATA_CSV.open(newline="", encoding="utf-8") as stream:
        return {row["species"]: (int(row["charge"]), int(row["multiplicity"]))
                for row in csv.DictReader(stream)}


# ---------- structure comparison ----------

def connectivity_mol(elements: list[str], coords: np.ndarray) -> Chem.Mol:
    """RDKit molecule with perceived connectivity only (no bond orders or charges)."""
    block = f"{len(elements)}\n\n" + "\n".join(
        f"{e} {x:.6f} {y:.6f} {z:.6f}" for e, (x, y, z) in zip(elements, coords))
    mol = Chem.MolFromXYZBlock(block)
    rdDetermineBonds.DetermineConnectivity(mol)
    return mol


def bonding_key(mol: Chem.Mol) -> str:
    """Canonical connectivity string including explicit H; equal keys = same bonding pattern."""
    editable = Chem.RWMol(mol)
    for atom in editable.GetAtoms():
        atom.SetNoImplicit(True)
        atom.SetNumRadicalElectrons(0)
    return Chem.MolToSmiles(editable, allHsExplicit=True, canonical=True)


def heavy_rmsd(probe: Chem.Mol, reference: Chem.Mol) -> float | None:
    """Heavy-atom RMSD after optimal alignment and atom mapping; None if the graphs don't match."""
    probe_heavy = Chem.RemoveHs(probe, sanitize=False)
    ref_heavy = Chem.RemoveHs(reference, sanitize=False)
    if ref_heavy.GetNumAtoms() < 2:
        return 0.0
    try:
        return rdMolAlign.GetBestRMS(probe_heavy, ref_heavy, maxMatches=MAX_MATCHES)
    except (RuntimeError, ValueError):
        return None


def n_linear_centres(elements: list[str], coords: np.ndarray) -> int:
    """Count two-coordinate heavy atoms whose bond angle is > LINEAR_ANGLE."""
    count = 0
    for centre, element in enumerate(elements):
        if element == "H":
            continue
        neighbours = [
            k for k in range(len(elements))
            if k != centre and np.linalg.norm(coords[k] - coords[centre])
            < 1.2 * (COVALENT_RADII[elements[k]] + COVALENT_RADII[element])
        ]
        if len(neighbours) != 2:
            continue
        u = coords[neighbours[0]] - coords[centre]
        v = coords[neighbours[1]] - coords[centre]
        cos_angle = np.clip(u @ v / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0)
        if np.degrees(np.arccos(cos_angle)) > LINEAR_ANGLE:
            count += 1
    return count


# ---------- analysis ----------

def categorize(best: dict) -> str:
    if not best["same_bonding"]:
        return "rearranged"
    if not best["converged"]:
        return "not_converged"
    if best["linear_artifact"]:
        return "linear_artifact"
    de = best["dE_vs_ref_kcal"]
    if de is None:
        return "no_reference_energy"
    if abs(de) <= DE_TOL:
        return "matches_reference"
    return "lower_than_reference" if de < 0 else "higher_than_reference"


def fmt(value, digits=4):
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return value


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: fmt(row.get(k)) for k in fields})


def analyze_model(model: str, metadata, references) -> tuple[list[dict], list[dict]]:
    model_dir = CONF_OPT_DIR / model
    conformers = []
    files = sorted(p for p in model_dir.glob("*_optimized.xyz") if CONF_FILE.match(p.name))
    for n, path in enumerate(files, 1):
        match = CONF_FILE.match(path.name)
        species, index = match["species"], int(match["index"])
        if species not in references:
            continue
        elements, coords, energy = read_xyz(path)
        steps, fmax_final = read_bfgs_log(model_dir / f"{species}_conf{index}.log")
        ref = references[species]
        mol = connectivity_mol(elements, coords)
        conformers.append({
            "model": model, "species": species, "conformer_index": index,
            "n_atoms": len(elements), "multiplicity": metadata[species][1],
            "energy_eV": energy, "n_steps": steps, "fmax_final": fmax_final,
            "converged": fmax_final is not None and fmax_final <= FMAX,
            "heavy_rmsd_vs_ref": heavy_rmsd(mol, ref["mol"]),
            "same_bonding": bonding_key(mol) == ref["bonding"],
            "linear_artifact": n_linear_centres(elements, coords) > ref["n_linear"],
        })
        if n % 100 == 0:
            print(f"  [{model}] {n}/{len(files)} conformers", flush=True)

    by_species = defaultdict(list)
    for row in conformers:
        by_species[row["species"]].append(row)

    best_rows = []
    for species, rows in sorted(by_species.items()):
        valid = [r for r in rows if r["energy_eV"] is not None]
        if not valid:
            continue
        best = min(valid, key=lambda r: r["energy_eV"])
        e_min = best["energy_eV"]
        for r in valid:
            r["rel_energy_kcal"] = (r["energy_eV"] - e_min) * EV_TO_KCAL

        ref_opt = REF_OPT_DIR / model / f"{species}_optimized.xyz"
        e_ref = read_xyz(ref_opt)[2] if ref_opt.exists() else None
        de = (e_min - e_ref) * EV_TO_KCAL if e_ref is not None else None
        for r in valid:
            r["dE_vs_ref_kcal"] = (r["energy_eV"] - e_ref) * EV_TO_KCAL if e_ref is not None else None

        rmsds = [r["heavy_rmsd_vs_ref"] for r in valid if r["heavy_rmsd_vs_ref"] is not None]
        summary = dict(best)
        summary.update({
            "n_conformers": len(valid),
            "best_conformer_index": best["conformer_index"],
            "dE_vs_ref_kcal": de,
            "min_heavy_rmsd_any_conformer": min(rmsds) if rmsds else None,
            "conformer_energy_span_kcal": max(r["rel_energy_kcal"] for r in valid),
        })
        summary["category"] = categorize(summary)
        best_rows.append(summary)
    return conformers, best_rows


def model_summary(model: str, conformers: list[dict], best_rows: list[dict]) -> dict:
    cats = Counter(r["category"] for r in best_rows)
    des = [r["dE_vs_ref_kcal"] for r in best_rows if r["dE_vs_ref_kcal"] is not None]
    rmsds = [r["heavy_rmsd_vs_ref"] for r in best_rows if r["heavy_rmsd_vs_ref"] is not None]
    n = len(best_rows)
    return {
        "model": model,
        "n_conformers": len(conformers),
        "n_conformers_converged": sum(r["converged"] for r in conformers),
        "n_species": n,
        "pct_matches_reference": 100.0 * cats["matches_reference"] / n if n else None,
        "median_dE_vs_ref_kcal": statistics.median(des) if des else None,
        "median_heavy_rmsd": statistics.median(rmsds) if rmsds else None,
        **{f"n_{c}": cats[c] for c in (
            "matches_reference", "lower_than_reference", "higher_than_reference",
            "linear_artifact", "not_converged", "rearranged", "no_reference_energy")},
    }


def species_consensus(all_best: list[dict]) -> list[dict]:
    by_species = defaultdict(list)
    for r in all_best:
        by_species[r["species"]].append(r)
    rows = []
    for species, rs in by_species.items():
        cats = Counter(r["category"] for r in rs)
        des = [r["dE_vs_ref_kcal"] for r in rs if r["dE_vs_ref_kcal"] is not None]
        rows.append({
            "species": species,
            "multiplicity": rs[0]["multiplicity"],
            "n_models": len(rs),
            "n_conformers": rs[0]["n_conformers"],
            "mean_dE_vs_ref_kcal": statistics.mean(des) if des else None,
            "min_dE_vs_ref_kcal": min(des) if des else None,
            "max_dE_vs_ref_kcal": max(des) if des else None,
            **{f"n_{c}": cats[c] for c in (
                "matches_reference", "lower_than_reference", "higher_than_reference",
                "linear_artifact", "not_converged", "rearranged")},
            "lower_in_models": ";".join(sorted(r["model"] for r in rs if r["category"] == "lower_than_reference")),
            "rearranged_in_models": ";".join(sorted(r["model"] for r in rs if r["category"] == "rearranged")),
        })
    rows.sort(key=lambda r: (int(r["species"].split("-")[0]), r["species"]))
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--models", nargs="+", default=None,
                        help="models to analyze (default: every folder in conformer_optimizations/)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    models = args.models or sorted(p.name for p in CONF_OPT_DIR.iterdir() if p.is_dir())
    metadata = load_metadata()

    references = {}
    for species in metadata:
        if species.endswith("TS"):
            continue
        elements, coords, _ = read_xyz(REFERENCE_DIR / f"{species}.xyz")
        mol = connectivity_mol(elements, coords)
        references[species] = {"mol": mol, "bonding": bonding_key(mol),
                               "n_linear": n_linear_centres(elements, coords)}

    OUTPUT_DIR.mkdir(exist_ok=True)
    all_conformers, all_best, summaries = [], [], []
    for model in models:
        if not (CONF_OPT_DIR / model).is_dir():
            print(f"Skipping {model}: no folder in {CONF_OPT_DIR}", file=sys.stderr)
            continue
        print(f"=== {model} ===", flush=True)
        conformers, best_rows = analyze_model(model, metadata, references)
        all_conformers += conformers
        all_best += best_rows
        summaries.append(model_summary(model, conformers, best_rows))

    conf_fields = ["model", "species", "conformer_index", "n_atoms", "multiplicity", "energy_eV",
                   "rel_energy_kcal", "dE_vs_ref_kcal", "n_steps", "fmax_final", "converged",
                   "heavy_rmsd_vs_ref", "same_bonding", "linear_artifact"]
    best_fields = ["model", "species", "multiplicity", "n_conformers", "best_conformer_index",
                   "dE_vs_ref_kcal", "heavy_rmsd_vs_ref", "min_heavy_rmsd_any_conformer",
                   "conformer_energy_span_kcal", "converged", "n_steps", "fmax_final",
                   "same_bonding", "linear_artifact", "category"]
    write_csv(OUTPUT_DIR / "conformer_table.csv", all_conformers, conf_fields)
    write_csv(OUTPUT_DIR / "best_conformers.csv", all_best, best_fields)
    write_csv(OUTPUT_DIR / "model_summary.csv", summaries, list(summaries[0]))
    consensus = species_consensus(all_best)
    write_csv(OUTPUT_DIR / "species_consensus.csv", consensus, list(consensus[0]))

    print(f"\n{'model':34s} {'% match':>8s} {'lower':>6s} {'higher':>7s} {'linear':>7s} {'noconv':>7s} {'rearr':>6s}")
    for s in summaries:
        print(f"{s['model']:34s} {s['pct_matches_reference']:8.1f} {s['n_lower_than_reference']:6d} "
              f"{s['n_higher_than_reference']:7d} {s['n_linear_artifact']:7d} {s['n_not_converged']:7d} "
              f"{s['n_rearranged']:6d}")
    print(f"\nWrote {OUTPUT_DIR}/conformer_table.csv, best_conformers.csv, model_summary.csv, species_consensus.csv")


if __name__ == "__main__":
    main()
