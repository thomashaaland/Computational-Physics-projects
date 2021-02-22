//   Project 3 a)
//   A first approach to calculating the 6D integral numerically
//   with Gauss-Legendre quadrature.

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gauleg.h"

//using namespace std;
std::ofstream ofile;


double gauss_legendre(double a, double b, double *r, double *u, int n) {
  double intgauss {0}; 
  gauleg::gauleg(a, b, r, u, n);
  
  //   evaluate the integral with the Gauss-Legendre method
  //   Note that we initialize the sum. Here brute force gauleg
  //   with same array, x, for all directions in 6D:
  
  for ( int i1 = 0;  i1 < n; i1++) {
    for ( int i2 = 0; i2 < n; i2++) {
      for (int i3 = 0; i3 < n; i3++) {
	for (int j1 = 0; j1 < n; j1++) {
	  for (int j2 = 0; j2 < n; j2++) {
	    for (int j3 = 0; j3 < n; j3++) {
	      double term = 1.0;
	      term *= gauleg::func_6D(r[i1],r[i2],r[i3],r[j1],r[j2],r[j3]);
	      term *= u[i1]*u[i2]*u[i3]*u[j1]*u[j2]*u[j3];
	      intgauss += term;
	    }
	  }
	}
      }
    }
  }
  return intgauss;
}


//   Main function begins here
int main()
{
  // Initial read in of some numbers:
  int n;
  int A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};//{10,15,20,25,30,35,40}; // Values of n to be tested 
  int size_of_A = sizeof(A)/sizeof(A[0]);
  
  double a, b; // Limits of the integral in Cartesian coordinates
  double const  pi = 3.14159265359;
  double exact_integral = 5*pi*pi/256; // Exact value of the integral to be evaluated
  double int_gauss = 0.0;
  double relatative_error = 0.0;
  double calculation_time = 0.0;
  std::string outfilename = "Gauss_Legendre_table.txt";
  //   reserve space in memory for vectors containing the mesh points
  //   weights and function values for the use of the gauss-legendre method
  double *r;
  double *u;

  std::cout << "Read in integration limits in the format: lower upper" << std::endl;
  std::cin >> a >> b;
  std::cout << "Got the limits " << a << " " << b << std::endl;

  std::cout << "Creating file " << outfilename << std::endl;
  ofile.open(outfilename);
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << " Integration limits in each dimension: a = " << a << " and b = " << b << std::endl;
  ofile << " n:        Result with Gauss-Legendre:      Exact result:       Relative error:     Calculation time [s]:" << std::endl;
  std::cout << "File created. Looping over " << size_of_A << " values" << std::endl;
  // loop over all values of n:
  for (int l=0; l < size_of_A; l++) {
    //inits
    std::cout << "\r" << (int)((double)(l+1)/(size_of_A)*100) << "% complete" << std::flush;
    n = A[l];
    r = new double [n];
    u = new double [n];
    int_gauss = 0.;
    calculation_time = 0.0;
    clock_t start, finish;
    start = clock();

    
    // trying the 6D integral with Gauss-Legendre method:
    int_gauss = gauss_legendre(a, b, r, u, n);

    finish = clock();
    calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;
    relatative_error = fabs(int_gauss - exact_integral)/exact_integral;
    
    ofile << std::setw(5) << std::setprecision(5) << n;
    ofile << std::setw(25) << std::setprecision(10) << int_gauss;
    ofile << std::setw(25) << std::setprecision(10) << exact_integral;
    ofile << std::setw(20) << std::setprecision(4) << relatative_error;
    ofile << std::setw(20) << std::setprecision(4) << calculation_time << std::endl;
    delete [] r;
    delete [] u;
  }
  std::cout << std::endl;

  ofile.close();

  return 0;
}  // end of main program


