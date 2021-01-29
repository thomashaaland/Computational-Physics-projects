# Jacobi method

The jacobi method diagonalises some matrix A such that A=R^T D R where D is diagonal. D in this way has the eigenvalues of matrix A, and R has the corresponding eigenvectors. To use this program you need to provide the matrix A and generate some matrix R. This program change the provided matrices inplace and returns void.

To run:
Run from main.cpp. The main file contain a small test to verify that the jacobi method diagonalises the A matrix. 
To use it in other programs, import jacobi.h and compile as normal, all functions in jacobi.cpp is under the namespace jacobi.
