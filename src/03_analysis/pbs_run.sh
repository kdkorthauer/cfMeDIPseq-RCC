#!/bin/bash
#PBS -l walltime=3:00:00,select=1:ncpus=1:mem=60gb
#PBS -N medip
#PBS -V
#PBS -A st-kdkortha-1
#PBS -o /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/JAN2020_medip_\$PBS_JOBID\.out
#PBS -e /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/JAN2020_medip_\$PBS_JOBID\.err

echo $summary

if [ $summary -eq 0 ]; then
  # Run this first
  R CMD BATCH --quiet --no-restore --no-save /arc/project/st-kdkortha-1/cfMeDIPseq-RCC/src/03_analysis/20200120_MEDIPS_RCC.R /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/20200120_MEDIPS_RCC.Rout
  R CMD BATCH --quiet --no-restore --no-save /arc/project/st-kdkortha-1/cfMeDIPseq-RCC/src/03_analysis/20200120_DMRs_RCC.R /scratch/st-kdkortha-1/_pbs/cfMeDIPseq//20200120_DMRs_RCC.Rout
else
  # Summarize
  R CMD BATCH --quiet --no-restore --no-save /arc/project/st-kdkortha-1/cfMeDIPseq-RCC/src/03_analysis20190730_RCC_figures_manuscript.R /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/20190730_RCC_figures_manuscript.Rout
fi