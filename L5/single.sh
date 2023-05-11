#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -t 15:00 --qos=short

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 2 ./E1

