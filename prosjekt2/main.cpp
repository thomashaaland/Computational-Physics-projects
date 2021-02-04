#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include "test.h"
#include "jacobi.h"

int main() {

    //Will solve the equation d^2/(dr^2) u = (r^2-L)u

    //define variables


    //declarations
    int k, l;
    k = 0;
    l = 0;
    int n = 150;
    int maxit = 100000;
    int iterations;
    double rmin = 0.0;
    double rmax = 4.37;
    double h = (rmax-rmin)/(n+1);
    double dh = h*h;
    double epsilon = 1e-8;
    arma::vec a = {0.01, 0.5, 1.0, 5.0};
    arma::vec r = arma::linspace(rmin, rmax, n+2);
    // R matrix holding all solutions
    arma::mat U = arma::zeros<arma::mat>(n,5);

    //constructing the tridiagonal matrix: d is diagonal element, e is off diagonal

    arma::vec d = arma::zeros(n);
    for(int i = 0 ; i<n ; i++) {
        d(i) = r(i+1)*r(i+1) + 2/dh;
    }
    arma::vec e = -1*arma::ones(n)/dh;

    //Create a matrix A, holding the info of the system, and the matrix B for finding solution
    arma::mat A = arma::zeros<arma::mat>(n,n);
    A.diag(+1)+=e(0);
    A.diag(-1)+=e(0);
    for(int i = 0 ; i<n ; i++) {
      A(i,i) = d(i);
    }
    int maxn = 50;
    int minn = 2;
    
    arma::mat to_file_iterations = arma::zeros(maxn-minn+1,2);
    for (int n = 2; n < maxn+1; n++){

        double h = (rmax-rmin)/(n+1);
        double dh = h*h;
	arma::vec d = arma::zeros(n);
	arma::vec r = arma::linspace(rmin, rmax, n+2);
        for(int i = 0 ; i<n ; i++) {
            d(i) = r(i+1)*r(i+1) + 2/dh;
        }

	arma::vec e = -1*arma::ones(n)/dh;
	arma::mat A = arma::zeros<arma::mat>(n,n);
        A.diag(+1)+=e(0);
        A.diag(-1)+=e(0);
        for(int i = 0 ; i<n ; i++) {
                A(i,i) = d(i);
        }

	arma::vec eigenvalues = jacobi::jacobi_eigenvalues(A, epsilon, maxit, k, l, n, iterations);
        to_file_iterations(n-minn,0) = n;
        to_file_iterations(n-minn, 1) = iterations-1;
	std::cout << n << std::endl;

    }
    to_file_iterations.save("iterations_n.dat", arma::raw_ascii);

    arma::vec eigenA;
    arma::mat eigvecA;
    eig_sym(eigenA, eigvecA, A);
    eigenA = sort(eigenA);
    for(int i = 0; i < n ; i++) {
        U(i,0) = eigenA(i);
    }

    arma::mat B = A;

    arma::vec eigenvalues = jacobi::jacobi_eigenvalues(B, epsilon, maxit, k, l, n, iterations);

    // For part c, create a new matrix A, containinga new potential V = ar^2+1/r


    arma::mat R1 = arma::zeros(n,n);
    arma::mat Rlist[] = {eigvecA, R1, R1, R1, R1};

    for(int k = 0 ; k < 4 ; k++) {
      arma::vec d2 = arma::zeros(n);
      for(int i = 0 ; i<n ; i++) {
	d2(i) = a(k)*a(k)*r(i+1)*r(i+1) + 1/r(i+1) + 2/dh;
      }
      
      for(int i = 0 ; i<n ; i++) {
	A(i,i) = d2(i);
      }
      eig_sym(eigenA,Rlist[k+1],A);
      eigenA = sort(eigenA);
      for(int j = 0; j < n ; j++) {
	U(j,k+1) = eigenA(j);
      }
    }
    
    arma::mat to_file1 = arma::zeros(n+2,5);
    for(int j=0 ; j < n ; j++) {
        to_file1(j,0) = U(j,0);
        to_file1(j,1) = U(j,1);
        to_file1(j,2) = U(j,2);
        to_file1(j,3) = U(j,3);
        to_file1(j,4) = U(j,4);
    }
    to_file1.save("solutions_u.dat", arma::raw_ascii);

    arma::mat to_file = arma::zeros(n+2,5);
    for(int j=0 ; j < n ; j++) {
        to_file(j,0) = Rlist[0](j,0);
        to_file(j,1) = Rlist[1](j,0);
        to_file(j,2) = Rlist[2](j,0);
        to_file(j,3) = Rlist[3](j,0);
        to_file(j,4) = Rlist[4](j,0);
    }
    to_file.save("eigvecs_u.dat", arma::raw_ascii);


    arma::mat C = B;
    arma::vec eigenvalues2 = jacobi::jacobi_eigenvalues(C, epsilon, maxit, k, l, n, iterations);
    std::cout << "C: " << eigenvalues2(0) << ",   " << eigenvalues2(1)
	      << ",   " << eigenvalues2(2) << std::endl;


    std::cout << "Victory!" << std::endl;
    test::UnitTests();
    return 0;

}


