# QM-GNNIS

## Publication

Katzberger P., Pultar F., and Riniker S., ChemRxiv. 2025

## Abstract

The conformational ensemble of a molecule is strongly influenced by the surrounding environment. Correctly modeling the effect of any given environment is, hence, of pivotal importance in computational studies. Machine learning (ML) has been shown to be able to model these interactions probabilistically, with successful applications demonstrated for classical molecular dynamics. While first instances of ML implicit solvents for quantum-mechanical (QM) calculations exist, the high computational cost of QM reference calculations hinders the development of a generally applicable ML implicit solvent model for QM calculations. Here, we present a novel way of developing such a general machine-learned QM implicit solvent model, which relies neither on QM reference calculations for training nor experimental data, by transferring knowledge obtained from classical interactions to QM. This strategy makes the obtained graph neural network (GNN) based implicit solvent model (termed QM-GNNIS) independent of the chosen functional and basis set. The performance of QM-GNNIS is validated on NMR and IR experiments, demonstrating that the approach can reproduce experimentally observed trends unattainable by state-of-the-art implicit-solvent models.

## Installation

````
# Install environment.
conda env create -f environment.yml
conda activate QMGNNIS

# Install repo after cloning from github
pip install .
````

To run the workflows you also need to install the [ORCA](https://www.faccts.de/orca/) quantum chemistry software.


## Usage

An example on how to use the tools are provided in the [demo.ipynb](notebooks/demo.ipynb) notebook.

## Reproducibility

This section is intended to provide a step-by-step guide to reproduce the results of the paper.

## Minimizations

### Molecular Balances
The minimizations for the 22 molecular balances were performed using the [run_minimisation_id_orca_min_hessian_no_mpi.py](production_runs/run_minimisation_id_orca_min_hessian_no_mpi.py) script. The submission script for each balance used for the minimizations are provided in the [submission_scripts](production_runs/submission_scripts/) folder. The results were analysed using the [analyse_Molecular_Balances_pub.ipynb](production_runs/analyse_Molecular_Balances_pub.ipynb) notebook.

### 2-methoxy-ethanol and 1,2-dimethoxyethane
The minimizations for the two compounds were performed using the [run_minimisation_id_orca_min_hessian_no_mpi_B3LYP_TZVP.py](production_runs/run_minimisation_id_orca_min_hessian_no_mpi_B3LYP_TZVP.py) script executed with the [submit_minimisations_B3LYP_TZVP.sh](production_runs/submit_minimisations_B3LYP_TZVP.sh) script. The results were analysed using the [analyse_Jtot_and_IR.ipynb](production_runs/analyse_Jtot_and_IR.ipynb) notebook.

## Authors
Paul Katzberger and Felix Pultar


