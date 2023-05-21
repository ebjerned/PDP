#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -n 16
#SBATCH -t 15:00 --qos=short

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 1 ./cg 1024
mpirun -np 4 ./cg 1024
mpirun -np 16 ./cg 1024

