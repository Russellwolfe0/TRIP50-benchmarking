# TRIP50 benchmark

[Website](https://russellwolfe0.github.io/TRIP50-benchmarking/) ·
[Final results](data/final/final_data_table.csv) ·
[Script reference](scripts/README.md)

TRIP50 compares DFT, machine-learned interatomic potentials (MLIPs), and g-xTB
on 50 triplet-state reactions. The current release contains 156 canonical
structures and 24 models: 7 DFT methods, 16 MLIPs, and g-xTB.

## Repository

```text
calculations/   Raw DFT outputs and MLIP/g-xTB result CSVs
configs/        Aliases, reaction types, charges, and multiplicities
data/           Intermediate and final benchmark tables
images/         Website and publication figures
scripts/        Calculation, parsing, analysis, and plotting tools
structures/     Canonical XYZ geometries
```

## Rebuild the benchmark

Python 3.10 or newer is recommended. The table pipeline uses only the standard
library; figures require Matplotlib:

```bash
python3 -m pip install -r requirements.txt
```

Regenerate tables from the committed calculation results:

```bash
# DFT
python3 scripts/dft_parser.py
python3 scripts/pivot_dft_results.py
python3 scripts/final_dft_results.py

# MLIPs
python3 scripts/mlip_parser.py
python3 scripts/pivot_mlip_results.py
python3 scripts/final_mlip_results.py

# g-xTB
python3 scripts/gxtb_parser.py
python3 scripts/pivot_gxtb_results.py
python3 scripts/final_gxtb_results.py

python3 scripts/final_data_table.py
python3 scripts/generate_plots.py
```

Each script overwrites only its documented derived outputs. Calculation drivers
are separate because rerunning them can be expensive.

## Run new calculations

Build charge and multiplicity metadata from the XYZ comments:

```bash
python3 scripts/chargemult_parser.py
```

Run an MLIP through the common ASE driver:

```bash
python3 scripts/run_mlip_single_points.py \
  --backend mace \
  --model MACE-OMOL \
  --output calculations/mlip_results/mace-omol.csv
```

Supported benchmark backends are AIMNet, FairChem, MACE, and ORB. Install ASE
and the selected model package in its own environment. Use `--dry-run` to verify
structure selection without loading the model.

Run g-xTB with the `gxtb` executable available on `PATH`:

```bash
python3 scripts/gxtb_driver.py --dry-run
python3 scripts/gxtb_driver.py
```

The driver gives each structure an isolated directory under
`calculations/gxtb/raw_outputs/` and writes its summary to
`calculations/gxtb_results/gxtb.csv`.

## Metrics

For total reactant energy `E_R`, transition-state energy `E_TS`, and total
product energy `E_P`:

```text
reaction energy = E_P  - E_R
forward barrier = E_TS - E_R
reverse barrier = E_TS - E_P
```

- Thermodynamic MAE is the mean absolute reaction-energy error.
- Per-reaction kinetic error is the mean absolute forward/reverse barrier error.
- Kinetic MAE averages that value over reactions.
- Combined MAE is `(thermodynamic MAE + kinetic MAE) / 2`.

Energies are reported in kcal mol⁻¹. Aliases in `configs/aliases.csv` always
resolve to their canonical structure before reaction quantities are calculated.

## Timing

`average_run_time_seconds` uses the most appropriate calculated time available
for each route:

| Route | Time used |
|---|---|
| Q-Chem | Reported CPU time |
| ORCA | Reported run time × parsed processor count |
| MLIP | CUDA-synchronized ASE energy-evaluation time |
| g-xTB | Timed `gxtb` subprocess |

Model loading is excluded. The values are not hardware-normalized, so comparisons
across CPU, GPU, and external-program routes should be treated as approximate.

## Main outputs

| Output | Contents |
|---|---|
| `data/final/final_data_table.csv` | Model category, runtime, and MAEs |
| `data/final/*_reaction_results.csv` | Reaction quantities and errors |
| `data/final/model_results_by_reaction_type.csv` | Per-model category results |
| `data/final/reaction_type_summary.csv` | Category-level summary |
| `data/final/*pareto*.csv` | Overall and category Pareto fronts |

See [scripts/README.md](scripts/README.md) for every script's inputs, outputs,
options, and pipeline position.
