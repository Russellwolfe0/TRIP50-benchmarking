# TRIP50 scripts

This directory contains the calculation drivers, parsers, aggregation tools, and
plot generators used by the TRIP50 benchmark. Unless noted otherwise, run every
command from the repository root:

```bash
python scripts/<script_name>.py
```

The scripts derive the repository root from their own location, so their default
file paths do not depend on the shell's current working directory. Running from
the root is still recommended because it makes custom relative paths predictable.

## Pipeline overview

The three calculation routes converge on the same reaction-level and model-level
tables:

```text
structures/reference/*.xyz
        |
        +--> chargemult_parser.py --> configs/structure_metadata.csv
        |
        +--> DFT output files
        |      dft_parser.py
        |        -> data/intermediate/dft_spe_results.csv
        |      pivot_dft_results.py
        |        -> data/intermediate/dft_species_pivot.csv
        |        -> data/final/dft_reaction_results.csv
        |      final_dft_results.py
        |        -> data/final/dft_benchmark_results.csv
        |
        +--> run_mlip_single_points.py
        |      -> calculations/mlip_results/*.csv
        |    mlip_parser.py
        |      -> data/intermediate/mlip_spe_results.csv
        |    pivot_mlip_results.py
        |      -> data/intermediate/mlip_species_pivot.csv
        |      -> data/final/mlip_reaction_results.csv
        |    final_mlip_results.py
        |      -> data/final/mlip_benchmark_results.csv
        |
        +--> gxtb_driver.py
               -> calculations/gxtb_results/gxtb.csv
             gxtb_parser.py
               -> data/intermediate/gxtb_spe_results.csv
             pivot_gxtb_results.py
               -> data/intermediate/gxtb_species_pivot.csv
               -> data/final/gxtb_reaction_results.csv
             final_gxtb_results.py
               -> data/final/gxtb_benchmark_results.csv

three benchmark tables
        |
        +--> final_data_table.py
               -> data/final/final_data_table.csv
        |
        +--> reaction_category_tables.py
               -> category CSVs and table images
        |
        +--> plot and website generators
```

The usual regeneration order is:

```bash
python scripts/chargemult_parser.py

python scripts/dft_parser.py
python scripts/pivot_dft_results.py
python scripts/final_dft_results.py

python scripts/mlip_parser.py
python scripts/pivot_mlip_results.py
python scripts/final_mlip_results.py

python scripts/gxtb_parser.py
python scripts/pivot_gxtb_results.py
python scripts/final_gxtb_results.py

python scripts/final_data_table.py
python scripts/reaction_category_tables.py
python scripts/pareto_front.py
python scripts/reaction_type_pareto_fronts.py
python scripts/reaction_type_pareto_fronts.py --absolute
python scripts/thermo_vs_kinetic_mae.py
python scripts/presentation_table.py
python scripts/generate_parity_plots.py
python scripts/generate_radar_plots.py
python scripts/generate_website_plots.py
```

Calculation drivers are intentionally omitted from that regeneration block:
rerunning them performs new model calculations and can be expensive.

## Timing conventions

The final model table deliberately uses calculated compute time rather than a
single undifferentiated wall-time field:

- **DFT:** `dft_parser.py` uses CPU time. Q-Chem reports CPU time directly.
  For ORCA, CPU time is estimated as total run time multiplied by the parsed
  `%pal nprocs` value. The provenance is retained in `cpu_time_source`.
- **MLIP:** `run_mlip_single_points.py` synchronizes CUDA immediately before and
  after `get_potential_energy()`. The measured interval is normalized by
  `mlip_parser.py` as `gpu_time_seconds` and propagated to the final table.
- **g-xTB:** `gxtb_driver.py` measures the `gxtb` subprocess calculation with
  `time.perf_counter()`. This becomes `calculation_time_seconds` downstream.

These values are useful for comparing methods within their routes. Cross-route
comparisons should still account for different hardware and accounting methods.

## Calculation and metadata scripts

### `chargemult_parser.py`

Builds the authoritative charge and spin-multiplicity table used by both
calculation drivers.

- **Inputs:** `structures/reference/*.xyz`, `configs/aliases.csv`, and
  `data/intermediate/reference_data.csv`.
- **Output:** `configs/structure_metadata.csv` with `species`, `charge`,
  `multiplicity`, `xyz_file`, and the original XYZ comment.
- **Behavior:** determines which species are physically required from the five
  reaction-state coefficients (`R1`, `R2`, `TS`, `P1`, `P2`), follows alias
  chains to canonical structures, and reads singlet/doublet/triplet labels from
  the second XYZ line. Every structure is assigned charge zero.
- **Validation:** rejects alias cycles, malformed XYZ files, missing or
  conflicting multiplicity labels, missing structures, and any canonical
  structure count other than 156.
- **CLI:** no arguments.

```bash
python scripts/chargemult_parser.py
```

### `run_mlip_single_points.py`

Universal ASE single-point driver for the supported MLIP backends.

- **Required options:** `--backend` and `--model`.
- **Optional inputs:** `--input` accepts one XYZ file or a directory and defaults
  to `structures/reference`; directories are searched recursively.
- **Output:** a standardized CSV. `--output` sets its location; otherwise the
  filename is derived from the model under `data/intermediate/mlip_results/`.
- **Metadata:** reads `configs/structure_metadata.csv`. `--charge` and
  `--multiplicity` can override the metadata for all selected structures.
- **Device:** `--device` defaults to `cuda`; CPU-capable backends may use
  `--device cpu`.
- **Timing:** one synchronized energy evaluation per structure is stored as
  `mean_time_seconds`. Model-loading time is excluded.
- **Failure handling:** writes a row for every attempted structure, records the
  exception in `error`, and returns exit status 1 if any calculation failed.
- **Dry run:** `--dry-run` lists selected structures and metadata without loading
  a model or writing results.

```bash
python scripts/run_mlip_single_points.py \
  --backend mace \
  --model MACE-OMOL \
  --output calculations/mlip_results/mace-omol.csv
```

Use `python scripts/run_mlip_single_points.py --help` for the current backend
choices and all options.

### `mlip_calculator_factories.py`

Internal adapter layer used by `run_mlip_single_points.py`; it is not normally
run directly.

- Lazily imports each MLIP package, which allows the common driver to exist in
  separate model-specific environments without importing every dependency.
- Maps the stable backend/model/device interface to ASE calculators.
- Supports shared calculators and per-structure calculators through
  `CalculatorProvider`.
- Contains adapters for FairChem, MACE, AIMNet, TorchANI, Orb, and MatterSim.
- Selects MACE POLAR versus OMOL from the model name and maps the public
  `MACE-OMOL` benchmark label to the checkpoint name expected by MACE.
- Raises targeted installation or model-name errors when a backend cannot be
  constructed.

`BACKENDS` is the authoritative source for the `--backend` choices exposed by
the MLIP driver.

### `gxtb_driver.py`

Runs organized g-xTB single-point calculations using the `gxtb` executable on
`PATH`.

- **Input:** `--input` accepts one XYZ file or a recursively searched directory;
  the default is `structures/reference`.
- **Metadata:** uses `configs/structure_metadata.csv` for charge and
  multiplicity.
- **Raw output:** `--work-dir` defaults to
  `calculations/gxtb/raw_outputs`. Each species receives its own directory with
  a copied XYZ, `.CHRG`, `.UHF`, and `gxtb.out`.
- **Summary output:** `--output` defaults to
  `calculations/gxtb_results/gxtb.csv` and follows the MLIP driver's general
  schema.
- **Timing:** measures the external calculation only and stores the interval in
  `mean_time_seconds`.
- **Options:** `--timeout SECONDS` applies per structure; `--dry-run` lists
  selected inputs without requiring the executable or writing results.
- **Failure handling:** preserves available output, records one failed CSV row,
  continues to the next structure, and exits 1 if any calculation failed.

```bash
python scripts/gxtb_driver.py
python scripts/gxtb_driver.py --input structures/reference/1-TS.xyz --dry-run
python scripts/gxtb_driver.py --timeout 600
```

## Parsing and normalization scripts

### `dft_parser.py`

Parses raw ORCA and Q-Chem output files into a common long-form species table.

- **Input:** all `*.out` files one directory below
  `calculations/single_points/raw_outputs/`; the directory name becomes the
  functional/model label.
- **Output:** `data/intermediate/dft_spe_results.csv`.
- **Species naming:** extracts a leading `<reaction>-R1`, `R2`, `TS`, `P1`, or
  `P2` token from each filename. Nonstandard names are skipped with a message.
- **ORCA parsing:** takes the last `FINAL SINGLE POINT ENERGY`, parses
  `TOTAL RUN TIME`, reads one-line or two-line `%pal nprocs`, and checks for
  normal termination.
- **Q-Chem parsing:** takes the last final-basis energy, reads reported wall and
  CPU times, and checks for the normal closing message.
- **Unknown formats:** retain a row but have blank numerical values and
  `completed=False`.
- **CLI:** no arguments.

### `mlip_parser.py`

Combines standardized MLIP result CSVs and converts them to the schema expected
by the pivot route.

- **Inputs:** positional CSV files and/or directories. Directories are searched
  recursively. With no arguments, it reads `calculations/mlip_results`.
- **Output:** `data/intermediate/mlip_spe_results.csv`, or the path supplied with
  `--output`.
- **Required input columns:** `species` and `energy_hartree`.
- **Model detection:** uses `model`, then `functional`, then the input filename.
- **Completion detection:** accepts common successful `status` or `completed`
  values; if neither exists, a nonblank energy is treated as complete.
- **Timing normalization:** chooses the first populated value among
  `gpu_time_seconds`, `mean_time_seconds`, `time_seconds`, and `time_s`, and
  writes it as `gpu_time_seconds`.
- **Safety:** excludes the output file itself from recursive input discovery.

```bash
python scripts/mlip_parser.py
python scripts/mlip_parser.py calculations/mlip_results/model.csv
python scripts/mlip_parser.py run_a/ run_b/ --output /tmp/mlip_spe_results.csv
```

### `gxtb_parser.py`

Normalizes one or more g-xTB driver CSVs for the shared pivot logic.

- **Default input:** `calculations/gxtb_results/gxtb.csv`.
- **Custom inputs:** positional CSV files or recursively searched directories.
- **Output:** `data/intermediate/gxtb_spe_results.csv`, configurable with
  `--output`.
- **Required columns:** `species`, `energy_hartree`, `mean_time_seconds`, and
  `status`.
- **Normalization:** writes `functional`, program `gxTB`, energy,
  `calculation_time_seconds`, completion state, and raw-output path.
- **Validation:** rejects missing columns, empty species, and duplicate
  model/species pairs.

```bash
python scripts/gxtb_parser.py
python scripts/gxtb_parser.py calculations/gxtb_results --output /tmp/gxtb.csv
```

## Pivot and reaction-energy scripts

### `pivot_dft_results.py`

Pivots completed DFT species calculations into reaction rows and calculates
thermodynamic and kinetic quantities relative to the coupled-cluster reference.

- **Inputs:** `data/intermediate/dft_spe_results.csv`, `configs/aliases.csv`, and
  `data/intermediate/reference_data.csv`.
- **Outputs:** `data/intermediate/dft_species_pivot.csv` and
  `data/final/dft_reaction_results.csv`.
- **Units:** converts Hartree to kcal/mol with 627.5095.
- **Aliases:** canonical calculations are authoritative; an alias receives the
  canonical species value even if an alias calculation exists.
- **Reaction definitions:** zero state coefficients in the reference table mean
  that a second reactant or product is absent, not that its energy is zero.
- **Quantities:** calculates reaction energy, forward barrier, reverse barrier,
  their signed errors, and kinetic error as the mean of the absolute forward
  and reverse barrier errors.
- **Timing:** carries species CPU time into state-specific pivot columns.
- **CLI:** no arguments.

### `pivot_mlip_results.py`

MLIP counterpart of `pivot_dft_results.py`.

- **Inputs:** `data/intermediate/mlip_spe_results.csv`, `configs/aliases.csv`, and
  `data/intermediate/reference_data.csv`.
- **Outputs:** `data/intermediate/mlip_species_pivot.csv` and
  `data/final/mlip_reaction_results.csv`.
- Uses the same alias handling, state-presence rules, unit conversion, reaction
  definitions, and error formulas as the DFT route.
- Carries `gpu_time_seconds` into state-specific `*_gpu_time_s` columns.
- Rejects duplicate completed calculations for a model/program/species key.
- **CLI:** no arguments.

### `pivot_gxtb_results.py`

Thin configuration wrapper around `pivot_mlip_results.py` for g-xTB.

- Changes the input to `data/intermediate/gxtb_spe_results.csv`.
- Writes `data/intermediate/gxtb_species_pivot.csv` and
  `data/final/gxtb_reaction_results.csv`.
- Maps `calculation_time_seconds` to state-specific
  `*_calculation_time_s` fields.
- All reaction calculations and validation come from the shared MLIP pivot
  implementation.
- **CLI:** no arguments.

## Benchmark aggregation scripts

### `final_dft_results.py`

Creates the compact DFT reaction benchmark table used by later summaries.

- **Inputs:** `data/intermediate/dft_species_pivot.csv` and
  `data/intermediate/reference_data.csv`.
- **Output:** `data/final/dft_benchmark_results.csv`.
- Recalculates thermochemistry and forward/reverse barriers for every reaction.
- Stores signed thermochemical and barrier errors.
- Defines per-reaction kinetic error as the mean absolute forward/reverse error.
- Defines `average_cpu_time_seconds` as the mean CPU time of the physically
  present reactants, products, and transition state.
- Missing required energies propagate as blank reaction quantities.
- **CLI:** no arguments; run `pivot_dft_results.py` first.

### `final_mlip_results.py`

Creates the corresponding compact MLIP reaction benchmark table.

- **Inputs:** `data/intermediate/mlip_species_pivot.csv` and the reference table.
- **Output:** `data/final/mlip_benchmark_results.csv`.
- Uses the same reaction and error formulas as `final_dft_results.py`.
- Averages the required structures' calculated GPU intervals into
  `average_gpu_time_seconds`.
- **CLI:** no arguments; run `pivot_mlip_results.py` first.

### `final_gxtb_results.py`

Thin wrapper around `final_mlip_results.py` for g-xTB.

- **Input:** `data/intermediate/gxtb_species_pivot.csv`.
- **Output:** `data/final/gxtb_benchmark_results.csv`.
- Reads state-specific calculation times and writes
  `average_calculation_time_seconds`.
- Uses the shared reaction and error implementation.
- **CLI:** no arguments; run `pivot_gxtb_results.py` first.

### `final_data_table.py`

Combines the DFT, MLIP, and g-xTB benchmark routes into the main model summary.

- **Inputs:** the three `*_benchmark_results.csv` files under `data/final`.
- **Output:** `data/final/final_data_table.csv`.
- Uses each route's intended time field: DFT CPU, MLIP GPU, and g-xTB
  calculation time.
- For each model, averages reaction-level run times and absolute thermochemical
  errors. Kinetic MAE is the mean of the existing per-reaction kinetic errors.
- Combined MAE gives thermochemical MAE and kinetic MAE equal weight:
  `(thermo_mae + kinetic_mae) / 2`.
- Rejects missing columns, models appearing in multiple method categories, and
  models with no values for a required aggregate.
- Sorts rows by method category and then model name.
- **CLI:** no arguments.

### `reaction_category_tables.py`

Calculates model performance separately for the seven reaction categories and
renders category tables.

- **Inputs:** `configs/rxntypes.csv`, the reference table, and the DFT, MLIP, and
  g-xTB species-pivot tables.
- **CSV outputs:** `data/final/model_results_by_reaction_type.csv` and
  `data/final/reaction_type_summary.csv`.
- **Image outputs:** per-category tables under `images/reaction_types/` and the
  category summary at `images/reaction_type_summary.png`.
- For each model/category pair, calculates average runtime, thermochemical MAE,
  kinetic MAE, and combined MAE from the category's reactions.
- The category summary reports category-level thermochemical MAE, kinetic MAE,
  and average time for DFT, MLIP, and semi-empirical models.
- Uses the same canonical aliases, present-state rules, and reaction formulas as
  the main pipeline.
- **CLI:** no arguments.

## Plot and website scripts

Plot scripts use Matplotlib's noninteractive `Agg` backend and can run on a
headless cluster. Install Matplotlib before running them.

### `pareto_front.py`

Creates the overall accuracy-versus-runtime Pareto plot.

- **Input:** `data/final/final_data_table.csv`, configurable with `--input`.
- **Outputs:** `images/pareto_front.png` and
  `data/final/pareto_front.csv`, configurable with `--output` and `--front-csv`.
- Minimizes positive average runtime and nonnegative combined MAE. A point is
  Pareto-optimal when no other model is at least as good on both axes and
  strictly better on one.
- Plots `log10(runtime)`, colors models by category, connects the front in
  increasing-runtime order, and labels only the models in `LABEL_MODELS`.
- `--dpi` controls output resolution and must be positive.

```bash
python scripts/pareto_front.py
python scripts/pareto_front.py --dpi 600 --output images/pareto_front_600dpi.png
```

### `reaction_type_pareto_fronts.py`

Creates separate thermodynamic and kinetic Pareto fronts for every reaction
category.

- **Input:** `data/final/model_results_by_reaction_type.csv`.
- **Default mode:** scales MAE within each category to a percentage of that
  category's model range and writes 14 plots to `images/Paretofrontsbycat/`
  plus `data/final/reaction_type_pareto_fronts.csv`.
- **Absolute mode:** `--absolute` retains kcal/mol MAE and writes a separate set
  under `images/Paretofrontsbycat_absolute/` plus
  `data/final/reaction_type_pareto_fronts_absolute.csv`.
- Runtime remains `log10(seconds)` in both modes.
- Pareto-front labels are vertically spaced to reduce overlap.
- **Current invariant:** the script requires exactly 24 model rows for each of
  the seven categories; update that validation if the benchmark model count is
  intentionally changed.

```bash
python scripts/reaction_type_pareto_fronts.py
python scripts/reaction_type_pareto_fronts.py --absolute
```

### `thermo_vs_kinetic_mae.py`

Plots model thermochemical MAE against kinetic MAE.

- **Input:** `data/final/final_data_table.csv`.
- **Output:** `images/thermo_vs_kinetic_mae.png` at 300 DPI.
- Uses equal axis limits and aspect ratio so the dashed `y = x` line is a true
  visual equality reference.
- Colors points by method category and labels only models listed in
  `LABEL_OFFSETS`.
- **CLI:** no arguments.

### `presentation_table.py`

Renders the complete final model table as a presentation-ready PNG.

- **Input:** `data/final/final_data_table.csv`.
- **Output:** `images/final_data_table.png` at 240 DPI.
- Sorts rows by combined MAE, formats runtimes as milliseconds, seconds,
  minutes, or hours, and rounds MAEs to two decimal places.
- Applies category colors and alternating row treatments; the image contains
  only the table.
- **Current invariant:** raises an error unless the input contains exactly 24
  models.
- **CLI:** no arguments.

### `generate_parity_plots.py`

Generates thermodynamic and kinetic reference-versus-prediction plots for every
model.

- **Default inputs:** the three `*_reaction_results.csv` files under
  `data/final`; positional arguments can replace them.
- **Output:** two files per model under `images/parity_plots/`, configurable with
  `--output-dir`.
- Thermodynamic plots contain one point per reaction. Kinetic plots contain both
  forward and reverse barriers with separate colors.
- Both use equal axes and a `y = x` line.
- Removes stale `*_parity.png` files from the selected output directory when
  their model is no longer present in the inputs.

```bash
python scripts/generate_parity_plots.py
python scripts/generate_parity_plots.py data/final/mlip_reaction_results.csv \
  --output-dir /tmp/parity
```

### `generate_radar_plots.py`

Creates one reaction-category radar chart per model for the website dropdown.

- **Input:** `data/final/model_results_by_reaction_type.csv`, configurable with
  `--input`.
- **Output:** `images/radar_plots/radar_<model-slug>.png`, configurable with
  `--output-dir`.
- Plots thermochemical and kinetic MAE across `C-C`, `C-O`, `C-S`, `HAT`,
  `Si-X`, `C-Hal`, and `N-X`.
- Requires exactly one value for every model/property/category combination.
- Automatically selects a radial maximum in 5 kcal/mol increments.
- Removes stale `radar_*.png` files not represented in the current input.

```bash
python scripts/generate_radar_plots.py
```

### `generate_website_plots.py`

Generates the additional fixed plots referenced by the archived-style website.

- **Scaled Pareto scatter:** reads `final_data_table.csv`, applies a square-root
  scaling to normalized combined MAE, and writes `images/paretov2_front.png`.
- **Method-type radar charts:** reads `reaction_type_summary.csv` and writes
  thermodynamic and kinetic DFT/MLIP/semi-empirical category comparisons under
  `images/radar_plots/`.
- **Reaction-coordinate plots:** reads all three reaction-result tables and
  writes `method_comparison.png` for the reactions in `SAMPLE_REACTIONS` under
  `images/reaction_coordinates/rxn<number>/`.
- The models used in reaction-coordinate comparisons are fixed in
  `SAMPLE_MODELS`. Change that constant if model names or the desired examples
  change.
- **CLI:** no arguments.

## Dependencies

All data-processing scripts use the Python standard library. Additional runtime
dependencies are needed only for calculations and graphics:

- **Plotting:** Matplotlib.
- **MLIP driver:** ASE plus the package for the selected backend; CUDA/PyTorch
  where required by that model.
- **g-xTB driver:** a working `gxtb` executable available on `PATH`.

Because model packages can have incompatible dependency requirements, the MLIP
factory imports packages lazily and is designed to be used from separate model
environments.

## Common maintenance points

- Model names propagate into CSVs, plot filenames, website dropdown values, and
  hard-coded label dictionaries. After renaming a model, regenerate downstream
  tables and plots and check `script.js`/`index.html` references.
- `configs/aliases.csv` is authoritative. Pivot scripts intentionally discard
  direct calculations for alias species in favor of the final canonical target.
- Do not interpret a zero species coefficient in `reference_data.csv` as a
  zero-energy calculation; it marks a species absent from that reaction.
- Several website plot generators contain fixed category, model, or reaction
  lists. Those constants are documented above and should be reviewed whenever
  benchmark membership changes.
- Parsers and aggregators overwrite their declared output files. Calculation
  drivers also overwrite their summary CSVs, while g-xTB species raw-output
  directories are reused.
