# TRIP50 benchmarking workflow

[View the TRIP50 benchmark website](https://russellwolfe0.github.io/TRIP50-benchmarking/)

This repository benchmarks electronic-structure methods on the TRIP50 reaction set. It provides a common workflow for density functional theory (DFT), machine-learned interatomic potentials (MLIPs), and g-xTB calculations, then evaluates their thermodynamic and kinetic errors producing tables and figures.

The assembled benchmark contains 50 reactions, 156 canonical molecular structures, and 24 models: 7 DFT methods, 16 MLIPs, and g-xTB. The main model-level result is [`data/final/final_data_table.csv`](data/final/final_data_table.csv), which reports each model's category, average calculation time, thermodynamic mean absolute error (MAE), kinetic MAE, and combined MAE.

## Repository layout

```text
TRIP50/
├── calculations/
│   ├── single_points/raw_outputs/  # Raw ORCA and Q-Chem DFT outputs
│   ├── mlip_results/               # CSV files written by the MLIP driver
│   ├── gxtb/raw_outputs/           # Organized g-xTB working/output directories
│   └── gxtb_results/               # g-xTB driver summary CSV
├── configs/
│   ├── aliases.csv                 # Original-to-canonical species mappings
│   ├── rxntypes.csv                # Reaction-category assignments
│   └── structure_metadata.csv      # Charge and multiplicity by structure
├── data/
│   ├── intermediate/               # Normalized single-point data and pivots
│   └── final/                      # Reaction- and model-level benchmark tables
├── images/                         # Generated tables and plots
├── scripts/                        # Drivers, parsers, aggregation, and plots
└── structures/reference/           # Reference XYZ geometries
```

All commands below assume the current directory is the repository root (`TRIP50/`).

## Requirements

The parsing and aggregation scripts use Python's standard library. Python 3.10 or newer is recommended. Plotting requires Matplotlib:

```bash
python3 -m pip install matplotlib
```

Running new MLIP calculations additionally requires ASE and the package associated with the selected backend. The driver supports AIMNet2, ANI, FairChem, MACE, MatterSim, and ORB calculators. Install only the backend required by the selected model, following that project's instructions.

Running g-xTB calculations requires a working `gxtb` executable. DFT jobs are not launched by this repository; completed ORCA or Q-Chem output files are parsed from the raw-output directory.

## Quick start: rebuild the existing results

Process the committed raw and driver results end to end with:

```bash
# DFT route
python3 scripts/dft_parser.py
python3 scripts/pivot_dft_results.py
python3 scripts/final_dft_results.py

# MLIP route
python3 scripts/mlip_parser.py
python3 scripts/pivot_mlip_results.py
python3 scripts/final_mlip_results.py

# g-xTB route
python3 scripts/gxtb_parser.py
python3 scripts/pivot_gxtb_results.py
python3 scripts/final_gxtb_results.py

# Combine all model categories
python3 scripts/final_data_table.py

# Category tables and figures
python3 scripts/reaction_category_tables.py
python3 scripts/pareto_front.py
python3 scripts/reaction_type_pareto_fronts.py
python3 scripts/reaction_type_pareto_fronts.py --absolute
python3 scripts/generate_radar_plots.py
python3 scripts/generate_parity_plots.py
python3 scripts/generate_website_plots.py
python3 scripts/presentation_table.py
python3 scripts/thermo_vs_kinetic_mae.py
```

These scripts overwrite their own derived CSV or PNG outputs, but do not alter raw calculation files.

## Structures and metadata

Canonical structure names follow `<reaction>-<state>`, where the reaction number is 1–50 and the state is:

| State | Meaning |
|---|---|
| `R1`, `R2` | First and second reactants |
| `TS` | Transition state |
| `P1`, `P2` | First and second products |

Some reactions have only one reactant or product. The reference table represents an absent second species with a zero contribution rather than a corresponding structure.

[`configs/aliases.csv`](configs/aliases.csv) maps historical structure names to canonical names. Every pivot route applies these aliases before matching structures to reactions. [`configs/structure_metadata.csv`](configs/structure_metadata.csv) is the authoritative source for neutral charge and spin multiplicity.

Regenerate metadata from XYZ comments and reference data with:

```bash
python3 scripts/chargemult_parser.py
```

The builder validates that all 156 canonical structures are present. The structure directory also contains `51-R1.xyz` and `51-TS.xyz`; drivers intentionally skip them because they are not part of the 50-reaction metadata.

## Calculation routes

### DFT

Place completed output files under a directory named for each functional:

```text
calculations/single_points/raw_outputs/
└── <functional>/
    ├── 1-R1.out
    ├── 1-TS.out
    └── ...
```

Run the route with:

```bash
python3 scripts/dft_parser.py
python3 scripts/pivot_dft_results.py
python3 scripts/final_dft_results.py
```

This produces, in order, normalized single-point data, species and reaction pivots, and the DFT benchmark table. The parser recognizes final energies and timing records from ORCA and Q-Chem outputs.

### MLIPs

[`scripts/run_mlip_single_points.py`](scripts/run_mlip_single_points.py) is a common ASE driver for every supported backend. It reads the canonical XYZ structures and charge/spin metadata, lazily loads the requested calculator, synchronizes CUDA around each calculation, and records energies and calculation times.

The required inputs are a backend and model identifier. For example:

```bash
python3 scripts/run_mlip_single_points.py \
  --backend mace \
  --model MACE-OMOL \
  --output calculations/mlip_results/mace-omol.csv \
  --device cuda
```

Verify structure discovery and metadata without loading a model:

```bash
python3 scripts/run_mlip_single_points.py \
  --backend mace \
  --model MACE-OMOL \
  --dry-run
```

Useful options are:

- `--input`: XYZ file or directory; default `structures/reference`
- `--output`: output CSV; place benchmark runs in `calculations/mlip_results/`
- `--device`: calculator device; default `cuda`
- `--dry-run`: validate inputs without running the calculator

Process the resulting CSV files with:

```bash
python3 scripts/mlip_parser.py
python3 scripts/pivot_mlip_results.py
python3 scripts/final_mlip_results.py
```

The parser discovers CSV files in `calculations/mlip_results/` by default. Files or directories can also be supplied as positional arguments.

### g-xTB

[`scripts/gxtb_driver.py`](scripts/gxtb_driver.py) runs each canonical XYZ file in a separate working directory. It copies the geometry, writes `.CHRG` and `.UHF`, captures program output, and creates a summary CSV consistent with the MLIP route.

First validate the planned calculations:

```bash
python3 scripts/gxtb_driver.py --dry-run
```

Then run it with `gxtb` available on `PATH`:

```bash
python3 scripts/gxtb_driver.py
```

The default input is `structures/reference`, metadata is read from `configs/structure_metadata.csv`, working directories are created under `calculations/gxtb/raw_outputs/`, and the summary is written to `calculations/gxtb_results/gxtb.csv`.

Process the summary with:

```bash
python3 scripts/gxtb_parser.py
python3 scripts/pivot_gxtb_results.py
python3 scripts/final_gxtb_results.py
```

## Benchmark definitions

For summed reactant energy `E_R`, transition-state energy `E_TS`, and summed product energy `E_P`:

```text
Thermodynamic reaction energy = E_P  - E_R
Forward activation energy     = E_TS - E_R
Reverse activation energy     = E_TS - E_P
```

Electronic energies are converted from hartree to kcal mol⁻¹ using 627.5095 where required.

- **Thermodynamic MAE** is the mean absolute reaction-energy error.
- **Kinetic error per reaction** is the mean of the absolute forward- and reverse-barrier errors.
- **Kinetic MAE** is the mean of those per-reaction kinetic errors.
- **Combined MAE** gives thermodynamic and kinetic performance equal weight: `(thermodynamic MAE + kinetic MAE) / 2`.

All final-table MAEs are in kcal mol⁻¹.

## Timing definitions and comparability

`average_run_time_seconds` contains calculated computation time rather than one uniform wall-clock metric, but its source differs by category:

| Category | Time used |
|---|---|
| DFT, Q-Chem | CPU time reported by Q-Chem |
| DFT, ORCA | Estimated CPU time: reported wall time × parsed processor count |
| MLIP | Mean synchronized calculation time around the ASE evaluation; model loading and warmup excluded |
| g-xTB | Mean elapsed subprocess calculation time |

These values support accuracy–cost visualization within this dataset, but they are **not hardware-normalized**. DFT CPU/core time, synchronized GPU evaluation time, and g-xTB process time are not equivalent. A rigorous speed comparison should rerun methods on specified hardware with consistent repetitions, resource accounting, precision, and inclusion rules.

## Final tables

| File | Contents |
|---|---|
| `data/final/dft_reaction_results.csv` | DFT quantities and errors by reaction and model |
| `data/final/mlip_reaction_results.csv` | MLIP quantities and errors by reaction and model |
| `data/final/gxtb_reaction_results.csv` | g-xTB quantities and errors by reaction |
| `data/final/*_benchmark_results.csv` | Route-specific reaction-level tables with timing |
| `data/final/final_data_table.csv` | Model category, time, and three MAE columns |
| `data/final/model_results_by_reaction_type.csv` | Model metrics split across seven categories |
| `data/final/reaction_type_summary.csv` | Category-level overall and model-class summaries |
| `data/final/pareto_front.csv` | Overall runtime/combined-MAE Pareto membership |
| `data/final/reaction_type_pareto_fronts.csv` | Scaled category Pareto membership |
| `data/final/reaction_type_pareto_fronts_absolute.csv` | Absolute-MAE category Pareto membership |

Reaction categories in [`configs/rxntypes.csv`](configs/rxntypes.csv) are C–C, C–O, C–S, hydrogen-atom transfer, Si–X, C–halogen, and N–X.

## Figures

| Command | Main output |
|---|---|
| `python3 scripts/presentation_table.py` | `images/final_data_table.png` |
| `python3 scripts/pareto_front.py` | `images/pareto_front.png` |
| `python3 scripts/thermo_vs_kinetic_mae.py` | `images/thermo_vs_kinetic_mae.png` |
| `python3 scripts/reaction_category_tables.py` | `images/reaction_type_summary.png` and `images/reaction_types/*.png` |
| `python3 scripts/reaction_type_pareto_fronts.py` | 14 scaled plots in `images/Paretofrontsbycat/` |
| `python3 scripts/reaction_type_pareto_fronts.py --absolute` | 14 absolute-MAE plots in `images/Paretofrontsbycat_absolute/` |
| `python3 scripts/generate_radar_plots.py` | One category-MAE radar plot per model in `images/radar_plots/` |
| `python3 scripts/generate_parity_plots.py` | Thermodynamic and kinetic parity plots in `images/parity_plots/` |
| `python3 scripts/generate_website_plots.py` | Scaled Pareto, method-type radar, and representative reaction-coordinate plots |

The overall Pareto plot uses average runtime on a logarithmic axis and combined MAE as its accuracy objective. Reaction-type Pareto plots are generated separately for thermodynamic and kinetic MAE.

## Validation and troubleshooting

- Start new MLIP and g-xTB campaigns with `--dry-run` to catch missing metadata or unexpected names.
- A failed calculation remains in the driver CSV with its status and error; inspect failed rows before pivoting.
- For missing species after normalization, check `configs/aliases.csv` and `configs/structure_metadata.csv` first.
- Pivot scripts require all species needed by a reaction and report incomplete model/reaction combinations.
- CUDA synchronization surrounds each measured MLIP calculation.
- Model identifiers are passed to backend libraries and may depend on the installed backend version.

## Reproducibility checklist

When reporting or extending these results, record:

1. model name, version, checkpoint, and numerical precision;
2. Python, ASE, backend-library, and g-xTB versions;
3. CPU/GPU model, processor count, and accelerator device;
4. whether initialization and I/O were included;
5. the calculation timing and hardware protocol;
6. the exact commit and changes to structures, aliases, metadata, or references.

# Citations

Found in the website documentation.
