ALFA-K takes longitudinal single cell sequencing data from an evolving
cell population as input, then estimates a local fitness landscape
encompassing thousands of karyotypes located near to the input data in
karyotype space. This repository contains source code and examples for
running ALFA-K, as well as an Agent Based Model (ABM) which simulates
evolving cell populations using fitness landscapes estimated by ALFA-K.

The repository is organized as follows:

The ABM requires compilation before use. E.g. to compile with the GCC
C++ compiler, change to ALFA-K root directory and run:

``` r
g++ ./ABM/main.cpp -o ./ABM/bin/ABM -std=c++17
```

For more details on the methods, see:

Beck, Richard J., and Noemi Andor. “Local Adaptive Mapping of Karyotype
Fitness Landscapes.” bioRxiv (2023): 2023-07.
