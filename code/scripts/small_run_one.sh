#!/usr/bin/R

module load R

#SBATCH --mem=100MB
#SBATCH --output=R-job.out
Rscript run_simulations.R 1 5 200 100 > out/small.one.out