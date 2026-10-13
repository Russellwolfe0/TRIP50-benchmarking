# WARNING: reference_smiles.csv was hand-corrected after this script made it
# (11-R1, 11-P2, 20-P1, 49-R1, 50-R2: azide/diazo/sulfonyl SMILES). Re-running
# this overwrites those fixes with the uncorrected structures/reference_smiles/*.smi.

from pathlib import Path
import csv
import os
import sys

REPO_DIR = Path(__file__).resolve().parents[1]  # scripts/ lives at the TRIP50 repo root
# Conformer-study data (raw runs and derived tables). Override with TRIP50_STUDY_DIR to point at
# a different copy, e.g. the large raw optimizations kept outside the repo.
STUDY_DIR = Path(os.environ.get("TRIP50_STUDY_DIR", REPO_DIR / "calculations" / "conformer_study"))
REFERENCE_SMILES_DIR = STUDY_DIR / "reference_smiles"
OUTPUT_CSV = STUDY_DIR / "reference_smiles.csv"

SPECIES_ORDER = {"R1": 0, "R2": 1, "P1": 2, "P2": 3}

def sort_key(row):
    rxn, role = row[0].split("-")
    return (int(rxn), SPECIES_ORDER[role])

if OUTPUT_CSV.exists() and "--force" not in sys.argv:
    sys.exit(f"{OUTPUT_CSV} already exists and is hand-corrected; pass --force to overwrite it.")

rows = []

for f in REFERENCE_SMILES_DIR.glob("*.smi"):
    species = f.stem
    if species.endswith("TS"):
        continue
    smiles = f.read_text().split()[0]

    rows.append((species, smiles))

rows.sort(key=sort_key)

with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(["code_name","SMILES"])
    writer.writerows(rows)
