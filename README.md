# Stress - network feedback in birds

[![DOI](https://zenodo.org/badge/193479297.svg)](https://zenodo.org/badge/latestdoi/193479297)


This repository contain all scripts and data needed to reproduce the analysis, figures and table presented in the manuscript: "Birds exposed to physiological stress post-breeding engage in stress-reducing social interactions in winter flocks." by Dekeukeleire et al

## How to

4 separate scripts are provided here:

- scripts/01\_create\_network\_indices.r: this script generate the winter network from the feeder data and derive the different network indices that will be used in the analysis.
- scripts/02\_analysis\_figure.r: this script runs the models and generate figure 3-5 from the manuscript
- scripts/03\_permutation\_analysis.r: this script runs the permutation analysis on the models
- scripts/04\_repeatability\_analysis.r: this script runs the repeatability analysis

The main data used here are the feeder recordings in data/feeder\_data\_formatted.csv and the stress and bird-level information (age, sex ...) in data/stress\_info.csv.

## Note on R and packages version 

All analysis presented in the manucript were ran with the following software and package versions:

- R v3.6.0
- rptR v0.9.22
- igraph v1.2.4.1
- asnipe v1.1.11

## License

The data posted here (in the data/ folder) are licensed under the [Creative Commons Attribution-NonCommercial 4.0 International license (CC BY-NC 4.0)](https://creativecommons.org/licenses/by-nc/4.0/), and all the R code is licensed under the [MIT license](LICENSE.md).
