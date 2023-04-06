for i in 1 2 4 5 8 10 16 20 25 32
do
	mpirun -np $i -bind-to none ./stencil ../../../../../home/maya/public/PDP_Assignment1/input1000000.txt test.txt 100 >> weak_scaling.txt
	echo "$i"
done
for i in 1 2 4 5 8 10 16 20 25 32
do
	mpirun -np $i -bind-to none ./stencil ../../../../../home/maya/public/PDP_Assignment1/input2000000.txt test.txt 100 >> weak_scaling.txt
	echo "$i"
done
for i in 1 2 4 5 8 10 16 20 25 32
do
	mpirun -np $i -bind-to none ./stencil ../../../../../home/maya/public/PDP_Assignment1/input4000000.txt test.txt 100 >> weak_scaling.txt
	echo "$i"
done
for i in 1 2 4 5 8 10 16 20 25 32
do
	mpirun -np $i -bind-to none ./stencil ../../../../../home/maya/public/PDP_Assignment1/input8000000.txt test.txt 100 >> weak_scaling.txt
	echo "$i"
done
