#!/bin/bash
#SBATCH --job-name=burnin
#SBATCH --output=burnin_%j.log
#SBATCH --error=burnin_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --mem=25G
#SBATCH --array=1-30

ml R
ml NLopt/2.7.0-GCCcore-11.2.0
ml LAPACK/3.12.1-GCC-11.3.0
ml GMP/6.2.1-GCCcore-12.2.0
ml MPFR/4.1.0-GCCcore-11.3.0
srun Rscript code/scenarios/burnin.R
