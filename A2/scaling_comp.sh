#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -n 25
#SBATCH -t 15:00 --qos=short
module unload intel
module unload intelmpi

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

for i in 16
do
    echo $i
    for j in input3600.txt input5716.txt input7488.txt input9072.txt
    do
        mpirun -np $i ./matmul /proj/uppmax2023-2-13/nobackup/matmul_indata/$j text.txt
    done
done

