# Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models

## Introduction

This repository contains the code for the paper ["Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models"](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1).

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