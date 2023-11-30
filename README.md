ALFA-K takes longitudinal single cell sequencing data from an evolving
cell population as input, then estimates a local fitness landscape
encompassing thousands of karyotypes located near to the input data in
karyotype space. This repository contains source code and examples for
running ALFA-K, as well as an Agent Based Model (ABM) which simulates
evolving cell populations using fitness landscapes estimated by ALFA-K.

The repository is organized as follows:

    ##                                                        levelName
    ## 1  ALFA-K                                                       
    ## 2   ¦--ABM                                                      
    ## 3   ¦   °--agent based model source code                        
    ## 4   ¦--examples                                                 
    ## 5   ¦   ¦--example_0                                            
    ## 6   ¦   ¦   °--visualizing karyotype frequencies and ALFA-K fits
    ## 7   ¦   ¦--example_1                                            
    ## 8   ¦   ¦   °--basics of running ABM and fitting landscapes     
    ## 9   ¦   °--example_2                                            
    ## 10  ¦       °--cross validation procedure for ALFA-K landscapes 
    ## 11  ¦--example_data                                             
    ## 12  ¦   °--example datasets                                     
    ## 13  °--utils                                                    
    ## 14      °--ALFA-K source code and utility functions

The ABM requires compilation before use. E.g. to compile with the GCC
C++ compiler, change to ALFA-K root directory and run:

``` r
g++ ./ABM/main.cpp -o ./ABM/bin/ABM -std=c++17
```

For more details on the methods, see:

Beck, Richard J., and Noemi Andor. “Local Adaptive Mapping of Karyotype
Fitness Landscapes.” bioRxiv (2023): 2023-07.
