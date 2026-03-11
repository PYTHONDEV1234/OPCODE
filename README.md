OPCODE

Oligodendrocyte Progenitor Cell Contrast-Oriented Detection Engine

Current Version: OPCCODE v1.0.0 (Research Freeze)
Status: Pre-publication validation complete

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

OPCODE is a contrast-based computational algorithm designed to identify oligodendrocyte progenitor cells (OPCs) in large-scale, single-cell RNA sequencing datasets of diverse mouse brain tissue.

The framework evaluates whether OPC-enriched genes are not only present in a cell but dominantly expressed as compared to competing genes associated with other cell types at the single-cell level to accurately identify whether cells are OPCs or not. This contrast-oriented scoring improves dynamic range and cluster identification of OPCs compared to conventional marker-expression averaging methods or module-scoring approaches which induce false positive detections and therefore create samples of OPCs containing higher impurity and less specificity relevant to demyelinating disease research.  

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Overview

Traditional OPC detection strategies include the averaging of expression rates of biologically proven genes associated with OPCs, or the application of gene module scoring through Scanpy. These methods confirm marker presence and can accurately identify OPC cells, but also have reduced cell-type specificity and samples of lower purity in that they do not assess whether OPC-associated genes are dominantly expressed over alternative cell-type programs, which better grades a cell as being an OPC or not. 

OPCODE introduces a contrast-based scoring system, where the following are calculated for every cell individually:

- Average expression of OPC-enriched genes

- Average expression of neuron-enriched genes

- Average expression of immune cell-enriched genes

The logic therefore flows like this:

- V5 Score = Mean(OPC genes) - Mean(Neuronal genes) - Mean(Immune genes).

This structure reduces false positives arising from cells showing weak marker expression of OPC genes which would otherwise be detected as OPC cells by current methods. 

V5 Score is measured without units, because it only calculates a relative transcriptional dominance index. Values of V5 below 0 indicate a dominant expression of non-OPC cell genes, around 0 indicate an inconclusive and likely non-OPC cell, and above 0 indicate a dominant expression of OPC-cell genes. The magnitude of V5 reflects the degree of expressional contrast within a dataset - between OPC and non-OPC associated genes - and is influenced by normalization strategy, expression variance, and biological composition rather than just by the size of the datasets. 

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Quick Start

1. Clone the repository

- git clone https://github.com/yourusername/opcode.git

2. Install required dependencies

- pip install -r requirements.txt

3. Run on a dataset

- python run_opc.py --input "path_to_your_dataset.h5ad" --output "outputs"

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Testing and Validation Datasets

OPCCODE has been evaluated on performance across numerous individual, independently-annotated mouse brain datasets:

- GSM1 (NCBI GEO)

- GSM2 (NCBI GEO)

- Macosko Mouse Brain Atlas

- 10Xv2 multi-region dataset (Allen Institute)

- 10Xv3 multi-region dataset (Allen Institute)

GSM1, GSM2, Macosko --> Independently clustered datasets used for validation and assessment of OPCODE's accuracy and precision. 

10Xv2, 10Xv3 --> Unlabeled, large-scale stress-test datasets (104 GB and 176 GB relatively, their scale makes them ideal for evaluating stability and scalability of OPCODE.)

Across all datasets, which differed in size, preprocessing and developmental methods, and labeling structure, OPCODE demonstrated:

- Increased dynamic range, relative to raw marker averaging and Scanpy module

- Improved cluster ranking accuracy

- Stable cross-region performance with various regions of the mouse brain

- Consistent high ranking of OPC clusters

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Validation Metrics

1. Dynamic Range:

- Difference between highest and lowest cluster mean V5 score. 

2. Cluster Rank Integrity Score (CRIS):

- Measures the normalized rank position of known OPC clusters among all clusters in a dataset, scaling between 0 (poor ranking) and 1 (optimal ranking). 

3. Cross-Region Stability:

- Consistency of OPC cluster ranking across various anatomical brain regions and datasets. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Limitations

- This algorithm was developed and validated only on public mouse brain datasets.

- The engine can only work on raw or processed .h5ad files for datasets.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Repository Structure

Root folder: opc_detection_tool/

- scoring engine/
	- v5_opc_scoring_engine.py

- analysis/
	- validation_metrics.py
	- method_comparison.py
	- full_benchmark_runner.py
	- stress_test_runner.py

- reporting/
	- pdf_report.py

- run_opc.py (Command-line interface)

- app.py (Streamlit interface)

- README.md (what you are reading right now)

----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Command-Line Usage

To run OPCODE from terminal:

python run_opc.py --input "path_to_your_dataset.h5ad" --output "outputs"

Outputs include:

- Scored dataset (.CSV)

- Summary statistics

- Optional validation metrics

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Requirements

- Python 3.10+

- scanpy

- pandas

- numpy

- matplotlib

- seaborn

- tqdm

To install dependencies:

pip install -r requirements.txt

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Versioning

v1.0.0 - Frozen during manuscript production and processing.

Algorithm and validation pipeline frozen for manuscript submission.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Data Availability

All datasets used for validation are publicly available

- GSM1 and GSM2: NCBI Gene Expression Omnibus

- Macosko Mouse Brain Atlas

- Allen Institute 10Xv2 and 10Xv3 Whole Brain Datasets

No other datasets were used in this study. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Code Availability

This repository is currently maintained in a pre-publication state

Full public release will occur upon publication of official manuscript

For academic inquiries or concerns, contact the author directly. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Reproducibility

OPCODE produces deterministic results after being applied to the same input dataset with normalized conditions. All validation metrics and benchmarking analyses can be reproduced using the scripts provided in the analysis directory and the instructions provided here to run the program. 

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

License

No open-source license has been applied at this stage.
All rights reserved.

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Author

Ansel Belani
anselbelani@gmail.com
March 2026