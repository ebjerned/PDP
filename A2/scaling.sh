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

for i in 1 4 9 16 25
do
    mpirun -np $i ./matrix /proj/uppmax2023-2-13/nobackup/matmul_indata/input10525.txt test.txt
done
for i in 1 4 9 16
do
    for j in input3600.txt input5716.txt input7488.txt input9072.txt input10525.txt
    do
        mpirun -np $i ./matrix /proj/uppmax2023-2-13/nobackup/matmul_indata/$j text.txt
    done
done

