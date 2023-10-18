#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=31200


module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

Rscript main.R $1 $2 $3 $4