#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=2g
#SBATCH -n 1
#SBATCH -t 96:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=<email>
module load r/4.1.3

#define variables
SUBDIR_NAME=N800/HCAR
NUM_PART=800
N_SIMS=5000

sbatch --output=/dev/null --error=/dev/null --time=10:00:00 --array=1-$N_SIMS --job-name=02_N_2000_Heap_Corr_HCAR --wait R CMD BATCH --no-save --no-restore "--args $SUBDIR_NAME $NUM_PART $N_SIMS" 02_Heaping_corr_sp_HCAR.R 02_Heaping_corr_sp_HCAR.Rout






