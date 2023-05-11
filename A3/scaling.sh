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


for k in 1 2 3
do
	for j in 2 32
	do
		for i in 125000000 250000000 500000000 1000000000
		do
			echo $k $j $i
			mpirun -np $j ./main /proj/uppmax2023-2-13/nobackup/qsort_indata/input$i.txt test.txt $k
		done
	done
done
