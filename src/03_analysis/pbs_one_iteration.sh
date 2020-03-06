#!/bin/bash
#PBS -l walltime=2:00:00,select=1:ncpus=1:mem=135gb
#PBS -N medip
#PBS -V
#PBS -A st-kdkortha-1
#PBS -o /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/JAN2020_medip_$PBS_JOBID\.out
#PBS -e /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/JAN2020_medip_$PBS_JOBID\.err
#PBS -M kdkorthauer@gmail.com

export WINDOWSIZE
export iter
export ntop

if [ $iter -eq 0 ]; then
  R CMD BATCH --quiet --no-restore --no-save /arc/project/st-kdkortha-1/cfMeDIPseq-RCC/src/03_analysis/20200120_MEDIPS_RCC.R /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/20200120_MEDIPS_RCC.Rout
else
  R CMD BATCH --quiet --no-restore --no-save /arc/project/st-kdkortha-1/cfMeDIPseq-RCC/src/03_analysis/20200120_DMRs_RCC.R /scratch/st-kdkortha-1/_pbs/cfMeDIPseq/20200120_DMRs_RCC_$iter\.Rout
fi
