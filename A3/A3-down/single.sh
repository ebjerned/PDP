#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 1
#SBATCH -A uppmax2023-2-36
#SBATCH -n 16
#SBATCH -t 30:00
module unload intel
module unload intelmpi

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 4 valgrind ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input93.txt test.txt 3

