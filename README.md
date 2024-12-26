# Rub_scVital - Single-Cell RNA Sequencing Integration Pipeline

## Overview

This Snakemake pipeline creates figures for scVital manuscript. The pipeline performs integration using our novel method, scVital, and comparison with five other integration methods. The goal is to provide a comprehensive analysis of cross-species cancer data integration.

## Features

- **Integration Methods**: Our novel integration method plus five other methods (scVital, Harmony, BBKNN, scVI, and scDREAMER).
- **Comparison Metrics**: Evaluation of integration quality using metrics such as ARI, FM, and LSS.
- **Visualization**: UMAP plots for integrated data visualization.
- **Reproducibility**: Fully reproducible analysis using Snakemake.

## Installation

1. **Clone the repository**:
    ```bash
    git clone https://github.com/betelab/Rub_scVital.git
    cd Rub_scVital
    ```

2. **Create and activate a conda environment**:
    ```bash
    conda create -n scInteg python=3.9
    conda activate scInteg
    ```

3. **Install dependencies**:
    ```bash
    conda install -c bioconda snakemake
    pip install scanpy anndata scikit-learn scVital
    ```

## Usage

1. **Configure the pipeline**:
    - Edit the `config/allParams.yaml` file to specify the input data and parameters for integration.

2. **Run the pipeline**:
    ```bash
    snakemake --cores 4 --use-conda
    ```

## Input Data

- **scRNA-seq Data**: Path to the scRNA-seq data file (e.g., `resources/data.h5ad`).

## Output

- **Integrated Data**: Integrated scRNA-seq data for mouse and human samples.
- **Comparison Results**: Metrics and visualizations comparing the integration methods.
- **Plots**: UMAP plots for visualizing the integrated data.

## Pipeline Structure

- `workflow/snakefile`: Defines the workflow of the pipeline.
- `resources/`: Directory where input files are stored.
- `config/allParams.yaml`: Configuration file for specifying input data and parameters.
- `worflow/scripts/`: Directory containing custom scripts for data processing and analysis.
- `worflow/notebooks/`: Directory containing custom ipthon notebook for figure generation.
- `results/`: Directory where output files will be saved.

## Contributing

We welcome contributions to improve this pipeline. Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

This README was made with the help of copilot.
