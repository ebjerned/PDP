#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -A uppmax2023-2-13
#SBATCH --reservation uppmax2023-2-13_3
#SBATCH -p core -n 10
#SBATCH -t 10:00


module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun ./E1 10000
