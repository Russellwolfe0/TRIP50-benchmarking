# Scripts

Run commands from the repository root. Defaults are repository-relative, and
scripts overwrite their declared outputs.

## Pipeline

| Stage | DFT | MLIP | g-xTB |
|---|---|---|---|
| Calculate | external | `run_mlip_single_points.py` | `gxtb_driver.py` |
| Parse | `dft_parser.py` | `mlip_parser.py` | `gxtb_parser.py` |
| Pivot | `pivot_dft_results.py` | `pivot_mlip_results.py` | `pivot_gxtb_results.py` |
| Aggregate | `final_dft_results.py` | `final_mlip_results.py` | `final_gxtb_results.py` |

After all three routes:

```bash
python3 scripts/final_data_table.py
python3 scripts/reaction_category_tables.py
python3 scripts/generate_plots.py
```

DFT aggregation uses CPU time, MLIP aggregation uses CUDA-synchronized
calculation time, and g-xTB uses its measured subprocess time.

## Calculation and metadata

### `chargemult_parser.py`

Reads `structures/reference`, `configs/aliases.csv`, and the reference data.
Writes `configs/structure_metadata.csv`. It resolves aliases, identifies required
reaction states, parses singlet/doublet/triplet labels from XYZ comments, and
requires all 156 neutral canonical structures. No options.

### `run_mlip_single_points.py`

Common ASE driver for AIMNet, FairChem, MACE, and ORB.

```bash
python3 scripts/run_mlip_single_points.py \
  --backend mace --model MACE-OMOL \
  --output calculations/mlip_results/mace-omol.csv
```

`--input` selects an XYZ or recursively searched directory, `--device` defaults
to CUDA, and `--charge`/`--multiplicity` override metadata. `--dry-run` validates
selection only. The output records energies, one synchronized evaluation time,
status, and errors; the command exits 1 if any calculation fails.

### `mlip_calculator_factories.py`

Internal lazy-loading adapters for the four supported MLIP packages. MACE-OMOL
and POLAR checkpoint selection and per-structure AIMNet calculators are handled
here. This module is imported by the MLIP driver, not run directly.

### `gxtb_driver.py`

Runs `gxtb --gxtb` for each selected XYZ using charge/multiplicity metadata.
Each species directory contains the XYZ, `.CHRG`, `.UHF`, and `gxtb.out`.

```bash
python3 scripts/gxtb_driver.py [--input PATH] [--work-dir PATH] \
  [--output FILE] [--timeout SECONDS] [--dry-run]
```

It continues after individual failures, records them in the summary CSV, and
exits 1 if the set is incomplete.

## Parsers

### `dft_parser.py`

Parses ORCA and Q-Chem `*.out` files beneath
`calculations/single_points/raw_outputs/<functional>/`. Writes energies,
completion status, timing, processor count, and provenance to
`data/intermediate/dft_spe_results.csv`. Q-Chem CPU time is reported directly;
ORCA CPU time is estimated from run time and `%pal nprocs`. No options.

### `mlip_parser.py`

Combines driver CSVs into `data/intermediate/mlip_spe_results.csv`. Positional
inputs may be files or recursively searched directories; the default is
`calculations/mlip_results`. `--output` changes the destination. It normalizes
model names, completion states, paths, and calculation time to
`gpu_time_seconds`.

### `gxtb_parser.py`

Normalizes g-xTB driver CSVs into `data/intermediate/gxtb_spe_results.csv`.
Accepts the same positional input and `--output` pattern as the MLIP parser.
Rejects missing fields and duplicate model/species rows.

## Reaction calculations

### `pivot_mlip_results.py`

Shared pivot implementation. It resolves aliases, converts Hartree to kcal/mol,
and calculates reaction energies, forward/reverse barriers, signed errors, and
mean absolute kinetic error. Incomplete reactions remain blank. Its MLIP defaults
write `mlip_species_pivot.csv` and `mlip_reaction_results.csv`.

### `pivot_dft_results.py` and `pivot_gxtb_results.py`

Thin wrappers around the shared pivot. They select the DFT CPU-time or g-xTB
calculation-time columns and route the corresponding input/output paths.

### `final_mlip_results.py`

Shared compact benchmark aggregation. It averages times across the physically
present states and writes reaction quantities and errors. Its defaults produce
`data/final/mlip_benchmark_results.csv`.

### `final_dft_results.py` and `final_gxtb_results.py`

Thin wrappers around the shared aggregation, selecting the route's time fields
and output file.

### `final_data_table.py`

Combines all three benchmark tables into `data/final/final_data_table.csv`.
For each model it reports average route time, thermodynamic MAE, kinetic MAE,
and their equally weighted combined MAE. No options.

### `reaction_category_tables.py`

Joins `configs/rxntypes.csv` to the pivot data and writes
`model_results_by_reaction_type.csv`, `reaction_type_summary.csv`, per-category
table images, and the summary image. No options.

## Figures

### `generate_plots.py`

The public plotting entry point. With no arguments it regenerates every table
image and plot. Pass one or more tasks to limit the work:

```bash
python3 scripts/generate_plots.py
python3 scripts/generate_plots.py pareto parity radar
```

Tasks are `tables`, `pareto`, `reaction-pareto`,
`reaction-pareto-absolute`, `accuracy`, `parity`, `radar`, and `website`.

The implementation modules are:

| Module | Output |
|---|---|
| `presentation_table.py` | Combined-MAE-sorted final table PNG |
| `pareto_front.py` | Runtime/combined-MAE front and CSV |
| `reaction_type_pareto_fronts.py` | Scaled or absolute category fronts |
| `thermo_vs_kinetic_mae.py` | Thermodynamic versus kinetic MAE |
| `generate_parity_plots.py` | Per-model thermodynamic/kinetic parity plots |
| `generate_radar_plots.py` | Per-model reaction-category radar plots |
| `generate_website_plots.py` | Scaled, method-type, and reaction-path plots |

These modules remain directly runnable when custom arguments are needed. Use
`python3 scripts/<name>.py --help` for their input/output options.

## Dependencies

Data processing uses the standard library. Figures require Matplotlib from the
root `requirements.txt`. Calculation drivers additionally require ASE, the
selected MLIP package, or the external `gxtb` executable.
