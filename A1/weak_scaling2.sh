for i in 1 2 4 5 8 10 16 20 25 32
do
	mpirun -np $i -bind-to none ./stencil ../../../../../home/maya/public/PDP_Assignment1/input100000000.txt test.txt 100 >> weak_scaling2.txt
	echo "$i"
done
