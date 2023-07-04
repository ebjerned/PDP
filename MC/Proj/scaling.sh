#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -n 32
#SBATCH -t 30:00 

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 1 ./mc 100000
mpirun -np 2 ./mc 100000
mpirun -np 4 ./mc 100000
mpirun -np 8 ./mc 100000
mpirun -np 16 ./mc 100000
mpirun -np 32 ./mc 100000

echo "large"
mpirun -np 16 ./mc 100000
mpirun -np 16 ./mc 160000
mpirun -np 16 ./mc 200000

echo "weak"
mpirun -np 1 ./mc 100000
mpirun -np 2 ./mc 200000
mpirun -np 4 ./mc 400000
mpirun -np 8 ./mc 800000
mpirun -np 16 ./mc 1600000

