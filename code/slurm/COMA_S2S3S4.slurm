#!/bin/bash
#SBATCH --job-name=COMA_S2S3S4
#SBATCH --output=COMA_S2S3S4_%j.log
#SBATCH --error=COMA_S2S3S4_%j.err
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --mem=150G
#SBATCH --array=1-30

ml R
ml NLopt/2.7.0-GCCcore-11.2.0
ml LAPACK/3.12.1-GCC-11.3.0
ml GMP/6.2.1-GCCcore-12.2.0
ml MPFR/4.1.0-GCCcore-11.3.0
srun Rscript code/scenarios/COMA_S2S3S4.R
