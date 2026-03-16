# OPCODE

**Oligodendrocyte Progenitor Cell Contrast-Oriented Detection Engine**

**Current Version:** OPCODE v1.1.0  
**Status:** Pre-publication validation complete

OPCODE is a contrast-based computational framework for identifying oligodendrocyte progenitor cells (OPCs) in large-scale single-cell RNA-sequencing datasets from diverse mouse brain tissue.

Rather than asking only whether OPC-associated genes are present, OPCODE evaluates whether the OPC transcriptional program is dominant over competing neuronal and immune programs at the single-cell level. This contrast-oriented strategy improves specificity, dynamic range, and biologically coherent cluster detection relative to conventional marker averaging or additive module scoring approaches, which can produce false positives and lower-purity candidate populations.

## Overview

Traditional OPC detection strategies often rely on:

- averaging expression of known OPC marker genes
- additive gene module scoring, for example Scanpy-style module scores

These approaches can detect OPC-associated signal, but they do not directly assess whether the OPC program is dominant relative to competing cell-type programs. As a result, they may identify cells with weak or mixed OPC-like signal as OPCs, lowering specificity.

OPCODE introduces a contrast-based scoring framework in which the following are calculated for each cell:

- average expression of OPC-enriched genes
- average expression of neuron-enriched genes
- average expression of immune cell-enriched genes

The core scoring logic is:

**V5 Score = Mean(OPC genes) - Mean(Neuronal genes) - Mean(Immune genes)**

This structure reduces false positives by penalizing cells in which non-OPC lineage programs remain strong even when some OPC markers are present.

## Interpreting the V5 Score

The V5 score is unitless. It functions as a relative transcriptional dominance index:

- **V5 < 0**: non-OPC programs dominate
- **V5 ~= 0**: mixed or inconclusive identity
- **V5 > 0**: OPC program dominates

The magnitude of V5 reflects transcriptional contrast within a dataset and can vary with:

- normalization strategy
- expression variance
- dataset composition
- biological heterogeneity

Accordingly, V5 is best interpreted primarily as a within-dataset ranking and separation metric, rather than a universal absolute scale across all datasets.

## Updates in v1.1.0

OPCODE v1.1.0 includes substantial workflow and robustness improvements over the earlier frozen baseline:

- normalization-aware preprocessing before V5 scoring
- improved cross-dataset score comparability
- automatic fallback generation of PCA, neighbors, and UMAP when missing
- automatic cluster fallback using Leiden when cluster metadata are absent
- improved PDF report layout and scientific wording
- improved handling of large datasets through sampled analysis mode
- improved gene-symbol normalization for Ensembl-style inputs
- improved progress feedback during analysis and report generation
- fixed experimental panel fallback text duplication
- improved visualization behavior for datasets with many clusters

## Quick Start

### 1. Clone the repository

    git clone https://github.com/PYTHONDEV1234/OPCODE.git
    cd OPCODE

### 2. Install dependencies

    pip install -r requirements.txt

### 3. Run the Streamlit app

    streamlit run app.py

### 4. Run from the command line

    python run_opc.py --input "path_to_your_dataset.h5ad" --output "outputs"

## Input Requirements

OPCODE currently supports:

- raw or processed `.h5ad` files
- mouse brain single-cell or single-nucleus RNA-seq datasets

Optional metadata can be merged when available.

If cluster labels are missing, OPCODE can automatically generate clustering structure from the input expression matrix.

## Validation Testing

OPCODE has been evaluated across multiple independently processed mouse brain datasets, including:

- GSE115746
- GSE60361
- GSM2906405 (Brain1)
- GSM2906406 (Brain2)
- Macosko Mouse Brain Atlas
- other larger datasets

Across datasets differing in:

- size
- preprocessing
- developmental composition
- labeling structure
- regional composition

OPCODE demonstrated:

- improved dynamic range relative to marker averaging and additive module scoring
- biologically coherent OPC cluster identification
- stable performance across mixed and OPC-enriched datasets
- the ability to detect both rare and abundant OPC populations
- repeated recovery of canonical oligodendrocyte-lineage markers in appropriate datasets

## Validation and Benchmarking Metrics

OPCODE includes several validation and benchmarking utilities.

- **Dynamic Range:** Difference between the highest and lowest cluster mean V5 score.
- **Cluster Rank Integrity Score (CRIS):** Measures how strongly the expected OPC population is prioritized relative to other clusters.
- **Cross-Region Stability:** Assesses consistency of OPC prioritization across anatomical regions and datasets.

Additional benchmarking scripts are included in the `analysis/` directory.

## Streamlit App Features

The Streamlit application (`app.py`) supports:

- upload of `.h5ad` datasets
- local file-path loading for large datasets
- automatic normalization when needed
- automatic UMAP generation when embeddings are missing
- automatic cluster generation when metadata are missing
- large-dataset sampled analysis mode
- scientific PDF report generation
- gene-symbol normalization for Ensembl-style gene identifiers

## Command-Line Usage

To run OPCODE from the terminal:

    python run_opc.py --input "path_to_your_dataset.h5ad" --output "outputs"

Typical outputs may include:

- scored dataset summaries
- validation metrics
- marker discovery results
- purification guidance
- scientific PDF report

## Repository Structure

    OPCODE/
    │   app.py
    │   README.md
    │   requirements.txt
    │   run_benchmark.py
    │   run_opc.py
    │
    ├── .streamlit/
    │   └── config.toml
    ├── analysis/
    ├── config/
    ├── reporting/
    ├── scoring_engine/
    └── utils/

## Key Components

- **app.py**: Streamlit interface
- **run_opc.py**: command-line interface
- **run_benchmark.py**: benchmark runner
- **analysis/**: benchmarking, validation, purification, visualization, and comparison tools
- **config/**: gene sets and surface-marker resources
- **reporting/**: PDF report generation and publication-oriented reporting tools
- **scoring_engine/**: V5 scoring implementation and core utilities
- **utils/**: metadata handling and helper utilities

## Requirements

Core dependencies include:

- Python 3.10+
- scanpy
- anndata
- pandas
- numpy
- scipy
- matplotlib
- seaborn
- tqdm
- reportlab
- streamlit

Install all dependencies with:

    pip install -r requirements.txt

## Limitations

- OPCODE has been developed and validated primarily on public mouse brain datasets
- the current workflow is centered on `.h5ad` input files
- interpretation of marker rankings should remain dataset-aware
- RNA enrichment does not guarantee protein-level sorting performance without experimental validation
- cross-dataset comparison of raw V5 magnitudes should be interpreted cautiously even after preprocessing standardization

## Reproducibility

OPCODE is designed to produce reproducible results on the same input dataset under the same preprocessing conditions.

The repository includes scripts for:

- validation
- benchmarking
- method comparison
- stress testing
- figure generation

These tools are intended to support reproducible manuscript preparation and downstream comparison studies.

## Code Availability

This repository reflects the pre-publication release state of OPCODE v1.1.0.

The codebase is publicly available here:

https://github.com/PYTHONDEV1234/OPCODE

## Data Availability

Validation datasets used in development and testing are publicly available from their original sources, including public GEO and atlas resources.

No proprietary datasets were used in the reported validation workflow.

## Versioning

- **v1.0.0**: frozen manuscript-era baseline
- **v1.0.1**: intermediate refinement stage
- **v1.1.0**: normalization-aware, report-improved, robustness-focused public release

## License

No open-source license has been applied at this stage.  
All rights reserved.

## Author

**Ansel Belani**

**Contact:**  
Institutional Email: abelani@mail.dccc.edu  
Personal Email: anselbelani@gmail.com  

March 2026