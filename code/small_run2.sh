#!/bin/bash
module load R/3.5.2
R CMD BATCH "--args 2 4 20 50" run_simulations.R out/tmp.try.2.out
