CONJUGATE GRADIENT METHOD
Author: Erik Bjerned, 2023-05-25
Parallel and Distributed Programming 1TD070


The conjugate gradient (CG) method is an algorithm used to solve linear systems of equations Au = b, 
under the prerequisite that A is symmetric and positive definite. The method is applied on a mesh of 
size n*n nodes. This program showcases a specific case of this method with predetermined A and b.

These instructions only apply to run this in general and will not consider systems as slurm.
The compilation and execution has been tested and verified on UPPMAX Snowy and Vitsippa.


REQUIERMENTS:
gcc
openmpi

USAGE:

./cg n

n - the number of nodes per side on a square mesh


OUTPUT:
The norm of the residual vector g.
(if #PRODUCE_OUTPUT is uncommented, the solution vector u will be 
printed to stdout after the final iteration. Rowmajor.)


