# Linear Programming for Protein Design

## Introduction

Welcome to the Linear Programming for Protein Design repository! This repository contains a Python package that designs diverse protein libraries by seeding linear programming with scores from inverse folding and protein language models. The software takes as input a matrix of *in silico* deep mutational scanning data, where each row corresponds to a mutation and each column corresponds to the score computed by a deep learning model. The software then uses linear programming to select a subset of mutations that maximize the diversity of the library while maintaining a high average score. The software outputs the selected mutations in a CSV file. The software is designed to be used by researchers in the field of protein engineering who want to design diverse protein libraries for directed evolution experiments. The software is open-source and can be freely downloaded and used by anyone.

This repository accompanies the paper ["Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models"](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1).

## Installation

Create an environment with Python >=3.9 and install the dependencies:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

- For developers:
```bash
black -S -t py39 lp_protein_design && isort lp_protein_design && flake8 --ignore=E501,E203,W503 lp_protein_design
```

## Run the code

```bash
lpsolve ./example_data/trastuzumab_spm.csv 10
```