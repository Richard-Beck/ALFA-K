# Chromosome Fitness Landscapes – alfak Analysis Code
ALFA-K takes longitudinal single cell sequencing data from an evolving
cell population as input, then estimates a local fitness landscape
encompassing thousands of karyotypes located near to the input data in
karyotype space. This repository contains source code and examples for
running ALFA-K, as well as an Agent Based Model (ABM) which simulates
evolving cell populations using fitness landscapes estimated by ALFA-K.
This repository contains all scripts and data needed to reproduce the analyses and figures from:

Beck, Richard J., and Noemi Andor. “Local Adaptive Mapping of Karyotype
Fitness Landscapes.” bioRxiv (2023): 2023-07. (pending citation update)

---

Quick Start
-----------

    ## clone the repo
    git clone https://github.com/Richard-Beck/ALFA-K.git
    cd ALFA-K
    ## install alfakR companion package:
    remotes::install_github("https://github.com/Richard-Beck/alfakR")
           
    # Run the pipeline
    Rscript scripts/01_preprocess_cell_lines.R
    # See note about batch scripts below
    ./scripts/02_batch_run_alfak.sh
    # ... etc.

All outputs will be written to `data/` and `figs/`.

---

Repository Structure
--------------------

.
├── README.md             # You're here  
├── scripts/              # Main analysis scripts  
├── R/                    # Utility functions (themes, packages, etc.)  
├── data/  
│   ├── raw/              # Immutable input data  
│   └── processed/        # Outputs from scripts   
├── figures/              # Generated figures 


---

Analysis Details
-----------------

Scripts should be run in order. 

Supplementary scripts (prepended with S) deal with generating and analyzing synthetic ABM data.

Some scripts are standalone, others should be run repeatedly with various data inputs. See the batch scripts for examples of how to set that up.

Running all the scripts as instructed will produce similar but not identical results to the companion article. In order to reproduce those results, extract the compressed raw data folder (/data/archive_raw.tar.gz) and ensure it is named 'raw'. This will supply the ABM synthetic simulations used in the manuscript. Download the processed data folder (22Gb), extract it, and ensure it is named 'processed'. This will provide the ALFA-K simulations used in the manuscript. Scripts 4-7, 9, and X01 should then reproduce the main figures - S04 and X02 produces the supplementary figures.






  



