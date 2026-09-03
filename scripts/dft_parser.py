'''Parses all of the files in TRIP50/calculations/single_points/raw_outputs. Assembles into 
CSV with columns "species", "functional", "program", "energy_hartree", "wall_time_seconds",
 "cpu_time_seconds", "cpu_time_source", "nprocs", "completed", "output_file"'''


from  pathlib import Path
import csv
import re

PROJECT_ROOT = Path(__file__).resolve().parents[1]

RAW_OUTPUTS = PROJECT_ROOT / "calculations/single_points/raw_outputs"
OUTPUT_CSV = PROJECT_ROOT / "data/intermediate/dft_spe_results.csv"

ORCA_ENERGY = re.compile(
    r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)"
)

ORCA_TIME = re.compile(
    r"TOTAL RUN TIME:\s*"
    r"(\d+)\s+days\s+"
    r"(\d+)\s+hours\s+"
    r"(\d+)\s+minutes\s+"
    r"(\d+)\s+seconds\s+"
    r"(\d+)\s+msec",
    re.IGNORECASE,
)

ORCA_NPROCS = re.compile(
    # ORCA may echo either ``%pal nprocs 8`` on one line or the two-line
    # form ``%pal`` followed by ``nprocs 8``. Support both.
    r"^\|\s*\d+>\s*%pal(?:\s+nprocs\s+|[\s\S]*?^\|\s*\d+>\s*nprocs\s+)(\d+)\b",
    re.IGNORECASE | re.MULTILINE,
)

QCHEM_ENERGY = re.compile(
    r"Total energy in the final basis set\s*=\s*(-?\d+\.\d+)"
)

QCHEM_TIME = re.compile(
    r"Total job time:\s*([\d.]+)s\(wall\),\s*([\d.]+)s\(cpu\)",
    re.IGNORECASE,
)

# The method portion can itself contain underscores (for example wB97M_V),
# so identify the standardized species at the beginning of the filename rather
# than splitting on the final underscore.
SPECIES_FROM_FILENAME = re.compile(r"^(\d+-(?:R[12]|P[12]|TS))(?:_|$)")

def identify_program(text):
    if "O   R   C   A" in text:
        return "ORCA"

    if "Welcome to Q-Chem" in text:
        return "Q-Chem"

    return "unknown"

def parse_orca(text):
    energy_matches = ORCA_ENERGY.findall(text)
    time_match = ORCA_TIME.search(text)

    nprocs_match = ORCA_NPROCS.search(text)
    nprocs = int(nprocs_match.group(1)) if nprocs_match else None

    energy = float(energy_matches[-1]) if energy_matches else None

    if time_match:
        days, hours, minutes, seconds, milliseconds = map(int, time_match.groups())
        wall_time = (
            days * 86400
            + hours * 3600
            + minutes * 60
            + seconds
            + milliseconds / 1000
        )
    else:
        wall_time = None

    cpu_time = wall_time * nprocs if wall_time is not None and nprocs is not None else None
    cpu_time_source = "estimated_from_wall_time_and_nprocs" if cpu_time is not None else None
    completed = "ORCA TERMINATED NORMALLY" in text

    return energy, wall_time, nprocs, cpu_time, cpu_time_source, completed


def parse_qchem(text):
    energy_matches = QCHEM_ENERGY.findall(text)
    time_match = QCHEM_TIME.search(text)
    energy = float(energy_matches[-1]) if energy_matches else None
    wall_time = float(time_match.group(1)) if time_match else None
    cpu_time = float(time_match.group(2)) if time_match else None
    nprocs = None
    cpu_time_source = "reported_by_qchem" if cpu_time is not None else None

    completed = "Thank you very much for using Q-Chem" in text

    return energy, wall_time, nprocs, cpu_time, cpu_time_source, completed

rows = []

for functional_folder in sorted(RAW_OUTPUTS.iterdir()):
    if not functional_folder.is_dir():
        continue

    functional = functional_folder.name

    for output_file in sorted(functional_folder.glob("*.out")):
        species_match = SPECIES_FROM_FILENAME.match(output_file.stem)
        if not species_match:
            print(f"Skipping nonstandard filename: {output_file.name}")
            continue

        species = species_match.group(1)

        text = output_file.read_text(encoding="utf-8", errors="replace")
        program = identify_program(text)

        if program == "ORCA":
            energy, wall_time, nprocs, cpu_time, cpu_time_source, completed = parse_orca(text)

        elif program == "Q-Chem":
            energy, wall_time, nprocs, cpu_time, cpu_time_source, completed = parse_qchem(text)

        else:
            energy = None
            wall_time = None
            nprocs = None
            cpu_time = None
            cpu_time_source = None
            completed = False

        rows.append(
            {
                "species": species,
                "functional": functional,
                "program": program,
                "energy_hartree": energy,
                "wall_time_seconds": wall_time,
                "cpu_time_seconds": cpu_time,
                "cpu_time_source": cpu_time_source,
                "nprocs": nprocs,
                "completed": completed,
                "output_file": str(output_file.relative_to(PROJECT_ROOT)),
            }
        )

OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as file:
    writer = csv.DictWriter(
        file,
        fieldnames=[
            "species",
            "functional",
            "program",
            "energy_hartree",
            "wall_time_seconds",
            "cpu_time_seconds",
            "cpu_time_source",
            "nprocs",
            "completed",
            "output_file",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {OUTPUT_CSV}")
