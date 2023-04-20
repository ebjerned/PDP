#!/bin/bash -l

#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13

#SBATCH -p core -n 2
#SBATCH -t 1:00

export OMPI_MCA_btl_openib_allow_ib=1
module load gcc openmpi
mpirun main
