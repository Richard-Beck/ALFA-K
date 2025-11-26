#! /bin/bash
#submit to Slurm
sbatch <<EOSUBMIT
#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --job-name=S01b_run_vary_misseg_abm_sims
#SBATCH --cpus-per-task=51
#SBATCH --mem-per-cpu=16G
#SBATCH --time=1-6
#SBATCH --qos=small
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err


# Load R module and check if Rscript is available
ml R/4.2.0
which Rscript  # Verify Rscript is available

Rscript scripts/S01b_vary_misseg_abm_sims.R

EOSUBMIT


