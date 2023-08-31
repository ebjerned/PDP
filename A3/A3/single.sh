#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 1
#SBATCH -A uppmax2023-2-13
#SBATCH -n 16
#SBATCH -t 15:00 --qos=short
module unload intel
module unload intelmpi

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 1


