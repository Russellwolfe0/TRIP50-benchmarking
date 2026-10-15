# TRIP50 benchmark

[Website](https://russellwolfe0.github.io/TRIP50-benchmarking/) ·
[Final results](data/final/final_data_table.csv) ·
[Script reference](scripts/README.md)

TRIP50 benchmarks DFT, machine-learned interatomic potentials (MLIPs) and g-xTB
on 50 triplet-state reactions against DLPNO-CCSD(T)/CBS. It covers 156 reference
structures and 24 methods: 7 DFT methods, 16 MLIPs and g-xTB. A second study
tests whether the MLIPs recover the reference conformers from RDKit starting
geometries.

## Repository

```text
calculations/
  single_points/      DFT outputs
  mlip_results/       MLIP single-point CSVs, one per model
  gxtb_results/       g-xTB single-point CSV
  conformer_study/    MLIP conformer optimizations, analysis and DFT checks
configs/              Aliases, reaction types, charges and multiplicities
data/intermediate/    Parsed and pivoted energies
data/final/           Benchmark tables
images/               Website and publication figures
scripts/              Drivers, parsers, analysis and plotting
structures/reference/ Reference XYZ geometries
index.html, script.js     Website: energy benchmark
conformers.html, conformers.js   Website: conformer study
style.css                 Shared website styles
```

## Energy benchmark

Install: `python3 -m pip install -r requirements.txt`. The table pipeline needs
only the standard library; figures need Matplotlib.

Rebuild everything from the committed results (run from the repository root):

```bash
for route in dft mlip gxtb; do
  python3 scripts/${route}_parser.py
  python3 scripts/pivot_${route}_results.py
  python3 scripts/final_${route}_results.py
done
python3 scripts/final_data_table.py
python3 scripts/reaction_category_tables.py
python3 scripts/generate_plots.py
```

New calculations (each in an environment with ASE and the model package):

```bash
python3 scripts/chargemult_parser.py   # charge/multiplicity from XYZ comments
python3 scripts/run_mlip_single_points.py --backend mace --model MACE-OMOL \
  --output calculations/mlip_results/mace-omol.csv   # AIMNet, FairChem, MACE, ORB
python3 scripts/gxtb_driver.py         # needs gxtb on PATH; --dry-run to check
```

**Metrics.** Reaction energy = E_P − E_R; forward barrier = E_TS − E_R; reverse
barrier = E_TS − E_P. Thermodynamic MAE is the mean absolute reaction-energy
error. Kinetic MAE is the mean over reactions of the mean absolute
forward/reverse barrier error. Combined MAE is their average. Energies are in
kcal mol⁻¹, and aliases in `configs/aliases.csv` resolve to their canonical
structure first.

**Timing.** Q-Chem: reported CPU time. ORCA: run time × processor count. MLIP:
CUDA-synchronized ASE evaluation time. g-xTB: timed subprocess. Model loading
is excluded, and times are not hardware-normalized.

**Main outputs** (`data/final/`): `final_data_table.csv` (runtime and MAEs per
method), `*_reaction_results.csv` (per-reaction errors),
`model_results_by_reaction_type.csv`, `reaction_type_summary.csv`, and the
`*pareto*.csv` fronts.

## Conformer study

For each of the 106 non-TS species, RDKit conformers (AQME CSEARCH, 553 in
total) are optimized with each of the 16 MLIPs. The lowest-energy conformer is
then compared with the same model's optimization of the reference geometry:
ΔE = E(best conformer) − E(optimized reference). A species matches when
|ΔE| ≤ 0.5 kcal mol⁻¹. MLIP energies are only compared as differences within
one model.

```text
calculations/conformer_study/
  reference_smiles.csv          SMILES per species (hand-corrected; see below)
  reference_smiles/             Open Babel SMILES from the reference XYZ
  reference_conformers/         AQME CSEARCH SDFs (log: CSEARCH_data.dat)
  conformer_xyz/                Multi-frame XYZ inputs to the optimizer
  conformer_optimizations/      Optimized conformers (.xyz) and BFGS logs per model
  optimizations/                Each model's optimization of the reference geometry
  conformer_optimization_results/  Per-model run CSVs
  conformer_analysis/           conformer_table, best_conformers, model_summary,
                                species_consensus
  best_conformers_ordered/      Best conformer per model, in reference atom order
  dft_checks/                   wB97M-V/def2-QZVPP single points: MLIP conformer
                                vs reference
```

Requires RDKit, NumPy and (for the optimizer) ASE:

```bash
python3 scripts/run_conf_mlip_optimizations.py --models uma-s-1p2   # GPU
python3 scripts/analyze_conformer_optimizations.py
python3 scripts/make_ordered_best_conformers.py
```

Set `TRIP50_STUDY_DIR` to use a copy of the study data outside the repository.
Trajectory files (`.traj`) are not committed.

`reference_smiles.csv` was corrected by hand after `make_reference_smiles_csv.py`
generated it. Open Babel assigns no formal charges, so azides (11-R1, 11-P2,
49-R1), the diazo ester 50-R2 and PhSO₂• (20-P1) needed charge-separated
SMILES. The script will not overwrite the file without `--force`. SMILES carry
no spin, so charge and multiplicity always come from
`configs/structure_metadata.csv`.

See [scripts/README.md](scripts/README.md) for each script's inputs, outputs
and options.
