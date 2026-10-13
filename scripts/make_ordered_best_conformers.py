#!/usr/bin/env python3
"""Write each model's lowest-energy conformer in the reference atom order.

For every model in conformer_optimizations/ and every species, picks the
lowest-energy optimized conformer (energy= tag in the extended XYZ, same rule
as analyze_conformer_optimizations.py), maps its atoms onto the DFT reference
by bonding graph (RDKit GetBestAlignmentTransform, symmetry-aware), and writes

    best_conformers_ordered/<model>/<species>_conf<i>_ordered.xyz

with the atoms in the same order as structures/reference/<species>.xyz and the
coordinates superimposed onto the reference, so the two overlay directly in
PyMOL (and pair_fit / rms_cur work by index). The comment line holds the model,
energy (eV) and heavy-atom / all-atom RMSD to the reference.

Species whose bonding differs from the reference (rearranged, or broken H
connectivity) can't be mapped and are listed at the end instead of written.

Run (in the TRIP50 or aqme env, which have RDKit):
    python make_ordered_best_conformers.py
"""

from __future__ import annotations

import re
from collections import defaultdict
import os
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolAlign

RDLogger.DisableLog("rdApp.*")

REPO_DIR = Path(__file__).resolve().parents[1]  # scripts/ lives at the TRIP50 repo root
# Conformer-study data (raw runs and derived tables). Override with TRIP50_STUDY_DIR to point at
# a different copy, e.g. the large raw optimizations kept outside the repo.
STUDY_DIR = Path(os.environ.get("TRIP50_STUDY_DIR", REPO_DIR / "calculations" / "conformer_study"))
CONF_OPT_DIR = STUDY_DIR / "conformer_optimizations"
REFERENCE_DIR = REPO_DIR / "structures" / "reference"
OUTPUT_DIR = STUDY_DIR / "best_conformers_ordered"

MAX_MATCHES = 20000  # symmetry-equivalent mappings tried
CONF_FILE = re.compile(r"^(?P<species>\d+-(?:R1|R2|P1|P2))_conf(?P<index>\d+)_optimized\.xyz$")
ENERGY = re.compile(r"(?<![\w])energy=(-?[\d.eE+-]+)")


def read_xyz(path: Path) -> tuple[Chem.Mol, float | None]:
    """RDKit mol (connectivity only) and energy in eV from an (extended) XYZ file."""
    lines = path.read_text().splitlines()
    n_atoms = int(lines[0])
    atoms = [" ".join(line.split()[:4]) for line in lines[2:2 + n_atoms]]
    mol = Chem.MolFromXYZBlock("\n".join([str(n_atoms), ""] + atoms))
    rdDetermineBonds.DetermineConnectivity(mol)
    match = ENERGY.search(lines[1])
    return mol, (float(match.group(1)) if match else None)


def main() -> None:
    references = {p.stem: read_xyz(p)[0] for p in REFERENCE_DIR.glob("*.xyz")}
    skipped = []
    n_written = 0

    for model_dir in sorted(d for d in CONF_OPT_DIR.iterdir() if d.is_dir()):
        model = model_dir.name
        candidates = defaultdict(list)  # species -> [(energy, index, path)]
        for path in model_dir.glob("*_optimized.xyz"):
            match = CONF_FILE.match(path.name)
            if match:
                candidates[match["species"]].append((int(match["index"]), path))

        out_dir = OUTPUT_DIR / model
        out_dir.mkdir(parents=True, exist_ok=True)

        for species in sorted(candidates):
            ref = references.get(species)
            if ref is None:
                skipped.append((model, species, "no reference xyz"))
                continue

            confs = []
            for index, path in candidates[species]:
                mol, energy = read_xyz(path)
                if energy is not None:
                    confs.append((energy, index, mol))
            if not confs:
                skipped.append((model, species, "no energies"))
                continue
            energy, index, mol = min(confs, key=lambda c: c[0])

            try:
                rms_all, _, atom_map = rdMolAlign.GetBestAlignmentTransform(
                    mol, ref, maxMatches=MAX_MATCHES)
                rms_heavy = rdMolAlign.GetBestRMS(
                    Chem.RemoveHs(mol, sanitize=False), Chem.RemoveHs(ref, sanitize=False),
                    maxMatches=MAX_MATCHES)
            except (RuntimeError, ValueError):
                skipped.append((model, species, f"conf{index}: bonding differs from reference"))
                continue

            rdMolAlign.AlignMol(mol, ref, atomMap=list(atom_map))  # superimpose onto reference
            positions = mol.GetConformer().GetPositions()
            probe_for_ref = {ref_i: probe_i for probe_i, ref_i in atom_map}

            lines = [str(ref.GetNumAtoms()),
                     f"{species} conf{index} {model} energy={energy:.6f} eV "
                     f"heavy_rmsd={rms_heavy:.4f} all_rmsd={rms_all:.4f}"]
            for ref_i in range(ref.GetNumAtoms()):
                probe_i = probe_for_ref[ref_i]
                x, y, z = positions[probe_i]
                symbol = mol.GetAtomWithIdx(probe_i).GetSymbol()
                lines.append(f"{symbol:<2s} {x:16.8f} {y:16.8f} {z:16.8f}")
            (out_dir / f"{species}_conf{index}_ordered.xyz").write_text("\n".join(lines) + "\n")
            n_written += 1

    print(f"Wrote {n_written} files to {OUTPUT_DIR}")
    if skipped:
        print(f"\nSkipped {len(skipped)}:")
        for model, species, reason in skipped:
            print(f"  {model:35s} {species:7s} {reason}")


if __name__ == "__main__":
    main()
