#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH --mem=1g
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mail-type=END,FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user= <email address>
module load r/4.1.3

#set number of sims once
N_SIMS=5000


SUBDIR_NAME=/Poisson/N800
simname=Poisson800
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=06_Compile_Sims_Poisson800 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 06_compile_sims.r 06_compile_sims_Poisson800.Rout

SUBDIR_NAME=/Nbin/N800
simname=Nbin800
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=06_Compile_Sims_Nbin800 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 06_compile_sims.r 06_compile_sims_Nbin800.Rout

SUBDIR_NAME=/ZIP/N800
simname=ZIP800
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=06_Compile_Sims_ZIP800 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 06_compile_sims.r 06_compile_sims_ZIP800.Rout

SUBDIR_NAME=/ZINB/N800
simname=ZINB800
sbatch --output=/dev/null --error=/dev/null --time=2:00:00 --job-name=06_Compile_Sims_ZINB800 --wait R CMD BATCH --no-save --no-restore "--args $N_SIMS $SUBDIR_NAME $simname" 06_compile_sims.r 06_compile_sims_ZINB800.Rout








