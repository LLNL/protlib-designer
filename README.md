# `protlib-designer` : Integer Linear Programming for Protein Library Design

![Status](https://img.shields.io/badge/Status-Active-green.svg)
![Python](https://img.shields.io/badge/Python-3.9-blue.svg)
[![Paper](https://img.shields.io/badge/Paper-Download-green.svg)](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1)

## Introduction

Welcome to the `protlib-designer` repository! This repository contains a Python package that designs diverse protein libraries by seeding linear programming with deep mutational scanning data (or any other data that can be represented as a matrix of scores per single-point mutation). The software takes as input the score matrix, where each row corresponds to a mutation and each column corresponds to a different source of scores, and outputs a subset of mutations that maximize the diversity of the library while Pareto-optimizing the scores from the different sources. 

The paper [Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1) uses this software to design diverse antibody libraries by seeding linear programming with scores computed by Protein Language Models (PLMs) and Inverse Folding models. 

<p align="center">
<img src="images/method_diagram.png" width="800">
</p>


## Getting Started

In this section, we provide instructions on how to install the software and run the code.

### Installation

Create an environment with Python >=3.9 and install the dependencies:
```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

If you want a nice development environment, you can install the development dependencies:
```bash
pip install -e .[dev]
```
which will allow you to run the tests and the linter. You can run the linting with:
```bash
black -S -t py39 protlib_designer && flake8 --ignore=E501,E203,W503 protlib_designer
```



### Run the code

To run the code to create a diverse protein library of size 10 from the example data, run the following command:

```bash
protlib-designer ./example_data/trastuzumab_spm.csv 10
```

We provide a rich set of command-line arguments to customize the behavior of `protlib-designer`. For example, the following command runs `protlib-designer` with a range of 3 to 5 mutations per sequence, enforcing the interleaving of the mutant order and balancing the mutant order, and using a weighted multi-objective optimization:

```bash
protlib-designer ./example_data/trastuzumab_spm.csv 10 --min-mut 3 --max-mut 5 --interleave-mutant-order True --force-mutant-order-balance True --weighted-multi-objective True
```


For more information on the command-line arguments, run:

```bash
protlib-designer --help
```

### Input data

The input to the software is a matrix of per-mutation scores (the csv file `trastuzumab_spm.csv` in the example above). Typically, the score matrix is defined by *in silico* deep mutational scanning data, where each row corresponds to a mutation and each column corresponds to the score computed by a deep learning model. See the example data in the `example_data` directory for an example of the input data format. The structure of the input data is shown below:

| MutationHL | score-1 | score-2 | ... | score-N |
|------------|--------|--------|-----|--------|
| AH106C     | -0.1    | 0.2    | ... | 0.3    |
| AH106D     | 0.2    | -0.3    | ... | -0.4    |
| ...        | ...    | ...    | ... | ...    |
| YH107A     | -0.3    | 0.4    | ... | -0.5    |
| ...        | ...    | ...    | ... | ...    |

Import notes about the input data:

1. The `MutationHL` column contains the mutation in the format : `WT_residue` + `chain` + `position_index` + `mutant_residue`. For example, `A+H+106+C = AH106C` represents the mutation of the residue at position 106 in chain H from alanine to cysteine.
2. The `score-1`, `score-2`, ..., `score-N` columns contain the scores computed by the deep learning models for each mutation. Typically, the scores are the negative log-likelihoods ratios of the mutant residue and the wild-type residue, computed by the deep learning model: 

    $$ s_{ij}^{\text{PLM}} =  -\log \left( \frac{p(x_i = a_j | w)}{p(x_i = w_i | w)} \right) =  -\log(p(x_i = a_j | w)) + \log(p(x_i = w_i | w)), $$

    where $w$ is the wild-type sequence, and $p(x_i = a_j | w)$ is the probability of the mutant residue $a_j$ at position $i$ given the wild-type sequence $w$ as estimated by a Protein Language Model (PLM) or an Inverse Folding model (or any other deep learning model). For example, in [Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models](https://www.biorxiv.org/content/10.1101/2024.11.03.621763v1), we used the scores computed by the [ProtBert](https://pubmed.ncbi.nlm.nih.gov/34232869/) and [AntiFold](https://arxiv.org/abs/2405.03370) models.

## Contributing

Please read [CONTRIBUTING.md](./CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Citing This Work

If you use this software in your research, please cite the following paper:

```latex
@article {Hayes2024.11.03.621763,
	author = {Hayes, Conor F. and Magana-Zook, Steven A. and Gon{\c c}alves, Andre and Solak, Ahmet Can and Faissol, Daniel and Landajuela, Mikel},
	title = {Antibody Library Design by Seeding Linear Programming with Inverse Folding and Protein Language Models},
	elocation-id = {2024.11.03.621763},
	year = {2024},
	doi = {10.1101/2024.11.03.621763},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2024/11/03/2024.11.03.621763},
	eprint = {https://www.biorxiv.org/content/early/2024/11/03/2024.11.03.621763.full.pdf},
	journal = {bioRxiv}
}
```

## License

`protlib-designer` is released under an MIT license. For more details, please see the
[LICENSE](./LICENSE) and [RELEASE](./RELEASE) files. All new contributions must be made under the MIT license.

SPDX-License-Identifier: MIT

LLNL-CODE-2001645