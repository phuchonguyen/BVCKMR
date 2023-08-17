#!/bin/bash
#SBATCH -J sim-bvckmr-03
#SBATCH -o sim-gp-output/slurm_%x_%a.out
#SBATCH -e sim-gp-output/slurm_%x_%a.err
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --array=1-100
#SBATCH --account=herringlab
#SBATCH -p herringlab,statdept-low,volfovskylab-low,common,scavenger
#SBATCH --mail-user=phn5@duke.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE

[ -d sim-gp-output/scenario-003 ] || mkdir sim-gp-output/scenario-003
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.0
module load R
R CMD BATCH BVCKMR_Sim03_RunCode.R