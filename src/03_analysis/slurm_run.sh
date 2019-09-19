#!/bin/bash
#SBATCH -J medips
#SBATCH -n 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 0-5:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH --mem=20G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../../_slurm/csaw_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../../_slurm/csaw_%j.err  # File to which STDERR will be written, %j inserts jobid

if [ $summary -eq 0 ]; then
  # Run this first
  R CMD BATCH --quiet --no-restore --no-save 20190712_MEDIPS_RCC.R log/20190712_MEDIPS_RCC.Rout
else
  # Summarize
  R CMD BATCH --quiet --no-restore --no-save 20190730_RCC_figures_manuscript.R log/20190730_RCC_figures_manuscript.Rout
fi