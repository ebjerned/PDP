mpirun -np 1 --bind-to none ./stencil ../../../../../../home/maya/public/PDP_Assignment1/input1000000.txt test.txt 100 >> weak.txt
mpirun -np 2 --bind-to none ./stencil ../../../../../../home/maya/public/PDP_Assignment1/input2000000.txt test.txt 100 >> weak.txt
mpirun -np 4 --bind-to none ./stencil ../../../../../../home/maya/public/PDP_Assignment1/input4000000.txt test.txt 100 >> weak.txt
mpirun -np 8 --bind-to none ./stencil ../../../../../../home/maya/public/PDP_Assignment1/input8000000.txt test.txt 100 >> weak.txt
