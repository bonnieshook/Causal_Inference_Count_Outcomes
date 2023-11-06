#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=<email address>
module load r/4.1.3

#define variables
SUBDIR_NAME=ZIP/N800/
NUM_PART=800
N_SIMS=5000

sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --array=1-$N_SIMS --job-name=04_N_800_ZIP --wait R CMD BATCH --no-save --no-restore "--args $SUBDIR_NAME $NUM_PART $N_SIMS" 04_ZIP.R 04_ZIP.Rout






