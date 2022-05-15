#!/bin/bash
#
#SBATCH --job-name=aall
#SBATCH --mail-user=xl110[AT]duke.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH -c 80
#SBATCH --output=slurm-R-job-%J.stdout
#SBATCH --error=slurm-R-job-%J.stderr

module load R/4.0.3-rhel8
R CMD BATCH AALL_small_Normcomb.R

