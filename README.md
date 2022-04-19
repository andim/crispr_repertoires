# A scaling law in CRISPR repertoire sizes arises from avoidance of autoimmunity

[![DOI](https://zenodo.org/badge/469239506.svg)](https://zenodo.org/badge/latestdoi/469239506)

This repository contains the source code associated with the manuscript

Hanrong Chen, Andreas Mayer, Vijay Balasubramanian: [A scaling law in CRISPR repertoire sizes arises from avoidance of autoimmunity](https://doi.org/10.1101/2021.01.04.425308), Current Biology 2022

It allows reproduction of the statistical analyses and numerical results reported in the manuscript.

## Installation requirements

Most code uses Julia v1.6.5. A number of standard scientific Julia packages are needed for the numerical simulations and visualizations: IJulia, NBInclude, PyPlot, SpecialFunctions, Distributions, DifferentialEquations, Plots, Clustering, StatsPlots. These can be installed using the Julia package manager.

The sequence bias analysis uses Python, and depends on the packages: numpy, scipy, matplotlib, pandas.

## Dataset download

The dataset used can be downloaded from the [CRISPRCasdb Download page](https://crisprcas.i2bc.paris-saclay.fr/Home/Download), specifically the "SQL dump" which is a PostgreSQL database backup. To restore this backup, first install PostgreSQL. Then create a database named "CRISPRCasdb" by running `createdb -U postgres CRISPRCasdb`, and restore the backup by running `psql -U postgres -d CRISPRCasdb < <dir>/20210121_ccpp.sql` on the command line. To create the tables used by the analysis codes, run the scripts in `scripts.sql` on pgAdmin.

## Contact

If you run into any difficulties running the code, please contact us at `andimscience@gmail.com` or `chen_hanrong@gis.a-star.edu.sg`.

## License

The source code is freely available under an MIT license. The plots are licensed under a Creative Commons attributions license (CC-BY).
