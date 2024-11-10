# Linear Programming for Protein Library Design

## Introduction

Welcome to the Linear Programming for Protein Library Design repository! This repository contains a Python package that designs diverse protein libraries by seeding linear programming with scores from inverse folding and protein language models. The software takes as input a matrix of *in silico* deep mutational scanning data, where each row corresponds to a mutation and each column corresponds to the score computed by a deep learning model. The software then uses linear programming to select a subset of mutations that maximize the diversity of the library while maintaining a high average score. The software outputs the selected mutations in a CSV file. The software is designed to be used by researchers in the field of protein engineering who want to design diverse protein libraries for directed evolution experiments. The software is open-source and can be freely downloaded and used by anyone.

This repository accompanies the paper [Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1).

<!-- add the image in /Users/landajuelala1/Code/abag/lp-protein-design/images/method_diagram.pdf -->
<p align="center">
<img src="images/method_diagram.png" width="800">
</p>


## Installation

Create an environment with Python >=3.9 and install the dependencies:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Run the code

To run the code to create a diverse protein library of size 10 from the example data, run the following command:

```bash
lpsolve ./example_data/trastuzumab_spm.csv 10
```

For more information on the command-line arguments, run:

```bash
lpsolve --help
```

## Input data

The key input to the software is a matrix of *in silico* deep mutational scanning data, where each row corresponds to a mutation and each column corresponds to the score computed by a deep learning model. See the example data in the `example_data` directory for an example of the input data format. The structure of the input data is shown below:

| MutationHL | score1 | score2 | ... | scoreN |
|------------|--------|--------|-----|--------|
| AH106C     | 0.1    | 0.2    | ... | 0.3    |
| AH106D     | 0.2    | 0.3    | ... | 0.4    |
| ...        | ...    | ...    | ... | ...    |
| YH107A     | 0.3    | 0.4    | ... | 0.5    |
| ...        | ...    | ...    | ... | ...    |


The `MutationHL` column contains the mutation in the format : `WT_residue` + `chain` + `position_index` + `mutant_residue`. For example, `A+H+106+C = AH106C` represents the mutation of the residue at position 106 in chain H from alanine to cysteine.

The `score1`, `score2`, ..., `scoreN` columns contain the scores computed by the deep learning models for each mutation. Typically, the scores are the negative log-likelihoods ratios of the mutant residue and the wild-type residue, computed by the deep learning model: 

$$ s_{ij}^{\text{PLM}} =  -\log \left( \frac{p(x_i = a_j | w)}{p(x_i = w_i | w)} \right) =  -\log(p(x_i = a_j | w)) + \log(p(x_i = w_i | w)), $$

where $w$ is the wild-type sequence, and $p(x_i = a_j | w)$ is the probability of the mutant residue $a_j$ at position $i$ given the wild-type sequence $w$ as estimated by a Protein Language Model (PLM) or an Inverse Folding model (or any other deep learning model). For example, in [Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1), we used the scores computed by the [ProtBert](https://pubmed.ncbi.nlm.nih.gov/34232869/) and [AntiFold](https://arxiv.org/abs/2405.03370) models.