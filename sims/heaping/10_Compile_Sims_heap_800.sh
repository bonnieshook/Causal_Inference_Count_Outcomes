#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=<email>
module load r/4.1.3

#set number of sims once
N_SIMS=5000


SUBDIR_NAME=/HeapCorr/N800/IH
simname=HeapCorr800_IH
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=03_Compile_Sims_heapcorr_800_IH --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 10_compile_sims_heap_exclude.r 10_compile_sims_heapcorr800_IH.Rout

SUBDIR_NAME=/HeapCorr/N800/HCAR
simname=HeapCorr800_HCAR
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=03_Compile_Sims_heapcorr_800_HCAR --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 10_compile_sims_heap_exclude.r 10_compile_sims_heapcorr800_HCAR.Rout

SUBDIR_NAME=/HeapInc/N800/IH
simname=HeapInc800_IH
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=03_Compile_Sims_heapinc_800_IH --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 10_compile_sims_heap_exclude.r 10_compile_sims_heapinc800_IH.Rout

SUBDIR_NAME=/HeapInc/N800/HCAR
simname=HeapInc800_HCAR
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=03_Compile_Sims_heapinc_800_HCAR --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 10_compile_sims_heap_exclude.r 10_compile_sims_heapinc800_HCAR.Rout


