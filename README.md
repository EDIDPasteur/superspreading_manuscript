# Quantifying superspreading in bacterial STI outbreaks using phylodynamics

This repository contains the data and scripts to reproduce the analyses and
figures of the manuscript:

> **Quantifying superspreading in bacterial STI outbreaks using phylodynamics.**


## Repository layout

```
.
├── analyses/
│   ├── hierarchical_model/   # Real-data (Neisseria) hierarchical BDMM-Prime analysis
│   └── simulations/
│       ├── trees/            # Simulate trees under a superspreading model (ReMaster → BDMM-Prime)
│       ├── sequences/        # Simulate alignments from the simulated trees (AliSim)
│       └── tree_statistics/  # Tree-shape summaries (Colless index, pairwise RF)
├── figures/
│   ├── paper_figures.R       # Reproduces every figure in the manuscript
│   └── data/                 # Pre-computed inputs to paper_figures.R
└── src/                      # Shared helpers (SLURM job runner, alnlen)
```

## System requirements

### Dependencies
- **BEAST 2** (2.7) with packages **BDMM-Prime** (≥ 2.4), **feast**, **ReMaster**, **SA**
- **R** (4.4.1) with:
  - `tidyverse` 2.0.0, `latex2exp` 0.9.8, `ggpubr` 0.6.0, `beastio` 0.3.3
  - `coda` 0.19-4.1, `ape` 5.8, `phangorn` 2.12.1, `phyloTop` 2.1.3
- **RStudio** 2024.04.2+764 (optional — scripts can also be run from the command line with `Rscript`)
- **Python** (≥ 3.10) with `biopython`, `dendropy`
- **IQ-TREE 3** (with AliSim)
- Read-processing toolchain (used in `analyses/hierarchical_model/process_reads.sh`):
  `sra-tools` 3.2.1, `fastqc` 0.12.1, `snippy` 4.6.0, `shovill` 1.1.0
- Phylogenetics toolchain (used in `analyses/hierarchical_model/process_clusters.sh`):
  `iqtree` 2.4.0, `gubbins` 3.4.3, `snp-sites` 2.5.1, `seqtk` 1.5
- A C compiler to build `src/alnlen.c`

### Tested on
macOS 26.3 with R 4.4.1 (RStudio 2024.04.2+764) and BEAST 2.7.

### Hardware
- **Reproducing the figures** from pre-computed data requires only a standard laptop (< 2 min).
- **Reproducing the full simulation study** requires a SLURM HPC cluster; `src/run_xml.sh` is a SLURM array-job script that submits hundreds of parallel BEAST runs.
- No non-standard hardware is required to run the helper scripts or build `alnlen`.

### Installation

```bash
# Build the alnlen helper (used by Simulate_alignments.sh)
gcc -O2 src/alnlen.c -o src/alnlen
```

Installing all required software from scratch (conda environments + BEAST packages) takes approximately 15–30 minutes on a standard internet connection.

## Reproducing the figures

The fastest path: regenerate every figure in the paper from the cached inputs in
`figures/data/`. To do that, open [`paper_figures.R`](figures/paper_figures.R) in RStudio and run it end-to-end. On a modern laptop this takes under 2 minutes.

The script produces five PDF files in `figures/`: `Figure1.pdf`, `Figure2.pdf`, `Figure3.pdf`, `Figure4.pdf`, and `Figure5.pdf`, corresponding directly to the figures in the manuscript.

Inputs consumed by `paper_figures.R`:

| File | Description |
| --- | --- |
| `Simulation_results_trees.csv` | Posterior summaries from BDMM-Prime fits to simulated trees |
| `Simulation_results_seqs.csv`  | Posterior summaries from BDMM-Prime fits to simulated alignments |
| `Simulation_results_no_ss.csv` | Posterior summaries from the no-superspreading control |
| `BDMMPrime_neisseria_superspreading.thinned2x.log.gz` | BEAST log file for the hierarchical _Neisseria gonorrhoeae_ analysis |
| `tip_dates.txt`, `last_sample_date.txt` | Sampling dates used for skyline plots |

## Reproducing the analyses

### 1. Simulation study

#### 1a. Simulate trees and fit BDMM-Prime

```bash
cd analyses/simulations/trees
python run_simulations.py \
    -rf ReMaster_master_script.xml \
    -bf BDMMprime_master_script.xml \
    -b /path/to/beast \
    -d /path/to/output_dir
# No-superspreading control:
python run_simulations_no_ss.py [no ss files]
```

`run_simulations.py` sweeps four (f, M) combinations — the superspreader
fraction `f ∈ {0.05, 0.25}` and the relative transmission multiplier
`M ∈ {2, 10}` — generating one ReMaster XML per cell and one BDMM-Prime
inference XML per simulated tree (with and without tip types).

Submit the resulting BDMM-Prime XMLs with `src/run_xml.sh` (SLURM array job;
auto-stops when min ESS ≥ 200), then aggregate results separately for trees and
sequences, passing the output CSV name as the first argument:

```bash
# Aggregate tree-inference logs
Rscript analyses/simulations/trees/gather_results.R \
    Simulation_results_trees.csv path/to/tree_logs/*.log

# Aggregate sequence-inference logs
Rscript analyses/simulations/trees/gather_results.R \
    Simulation_results_seqs.csv path/to/seq_logs/*.log
```

#### 1b. Simulate alignments from the simulated trees

```bash
cd analyses/simulations/sequences
./Simulate_alignments.sh <xml_folder> <output_dir> ref.fa
```

Trees are extracted from the xml files and rescaled to substitutions/site by `rescale_tree.py` (rate
`4.5e-6 sub/site/year`), then AliSim simulates alignments under
`HKY{6.721}+I{0.9991}+G{0.5002}` rooted on `ref.fa`. The resulting alignments
are then fed back through BDMM-Prime as in 1a.


### 2. Hierarchical _Neisseria gonorrhoeae_ analysis

#### 2a. From reads to per-cluster alignments

```bash
cd analyses/hierarchical_model

# Per accession: download, QC, map, assemble
./process_reads.sh -a SRR [SRR...] -r reads/ -m mapping/ -s assembly/ -p reports/ -f reference.fa

# Per cluster: subset alignment, recombination filtering, ML tree
./process_clusters.sh -c cluster.clust -f full_alignment.fasta -o clusters_out/

# Append decimal sampling dates to FASTA headers
python rename_alignments.py -f cluster.fasta -o cluster.dated.fasta -d sample_dates.pkl
```

The pre-built per-cluster alignments and trees used in the paper are
provided as `alignments_2.tar.gz` and `trees.tar.gz`.

#### 2b. Hierarchical BDMM-Prime inference

```bash
beast -threads 2 -overwrite BDMMPrime_neisseria_hierarchical.xml
```

Alternative parametrizations explored live in
`different_parametrizations/` (`Default_param.xml`, `Default_sampprop.xml`,
`semi_reg.xml`, `semireg_multiplier.xml`, `full_reg.xml`).

