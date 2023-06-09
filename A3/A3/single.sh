#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -n 32
#SBATCH -t 15:00 --qos=short
module unload intel
module unload intelmpi

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

mpirun -np 2 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 3
mpirun -np 4 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 3
mpirun -np 8 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 3
mpirun -np 16 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 3
mpirun -np 32 ./quicksort /proj/uppmax2023-2-13/nobackup/qsort_indata/input125000000.txt test.txt 3


