#!/bin/bash -l
#SBATCH -M snowy
#SBATCH -p node -n 2
#SBATCH -A uppmax2023-2-13
#SBATCH -n 25
#SBATCH -t 15:00 --qos=short

module load gcc openmpi
export OMPI_MCA_btl_openib_allow_ib=1

echo 1200
mpirun -np 1 ./cg 1200
mpirun -np 4 ./cg 1200
mpirun -np 9 ./cg 1200
mpirun -np 16 ./cg 1200
mpirun -np 25 ./cg 1200

echo 2400
mpirun -np 1 ./cg 2400
mpirun -np 4 ./cg 2400
mpirun -np 9 ./cg 2400
mpirun -np 16 ./cg 2400
mpirun -np 25 ./cg 2400

echo 4800
mpirun -np 1 ./cg 4800
mpirun -np 4 ./cg 4800
mpirun -np 9 ./cg 4800
mpirun -np 16 ./cg 4800
mpirun -np 25 ./cg 4800

echo 512
mpirun -np 1 ./cg 512
mpirun -np 4 ./cg 1024
mpirun -np 9 ./cg 1536
mpirun -np 16 ./cg 2048
mpirun -np 25 ./cg 2560

echo 1024
mpirun -np 1 ./cg 1024
mpirun -np 4 ./cg 2048
mpirun -np 9 ./cg 3072
mpirun -np 16 ./cg 4096
mpirun -np 25 ./cg 5120
