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
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define epsilon 1.0E-6

//using namespace std;
std::ofstream ofile;

//     Here we define various functions called by the main program:
void gauss_laguerre(double *, double *, int, double);
void gauleg(double, double, double *, double *, int);
double func_6D(double, double, double, double, double, double);
double abs_dist(double, double, double, double, double, double);

//   Main function begins here
int main()
{
  // Initial read in of some numbers:
  int n;
  int A[] = {10,15,20,25,30,35,40}; // Values of n to be tested 
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

  std::cout << "Read in integration limits" << std::endl;
  std::cin >> a >> b;

  ofile.open(outfilename);
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << " Integration limits in each dimension: a = " << a << " and b = " << b << std::endl;
  ofile << " n:        Result with Gauss-Legendre:      Exact result:       Relative error:     Calculation time [s]:" << std::endl;
  // loop over all values of n:
  for (int l=0; l < 7; l++) {

  n = A[l];
  // trying the 6D integral with Gauss-Legendre method:
  r = new double [n];
  u = new double [n];
  int_gauss = 0.;
  calculation_time = 0.0;

  gauleg(a, b, r, u, n);

  //   evaluate the integral with the Gauss-Legendre method
  //   Note that we initialize the sum. Here brute force gauleg
  //   with same array, x, for all directions in 6D:
  clock_t start, finish;
  start = clock();
  for ( int i1 = 0;  i1 < n; i1++) {
    for ( int i2 = 0; i2 < n; i2++) {
      for (int i3 = 0; i3 < n; i3++) {
        for (int j1 = 0; j1 < n; j1++) {
          for (int j2 = 0; j2 < n; j2++) {
            for (int j3 = 0; j3 < n; j3++) {
              double term = 1.0;
              term *= func_6D(r[i1],r[i2],r[i3],r[j1],r[j2],r[j3]);
              term *= u[i1]*u[i2]*u[i3]*u[j1]*u[j2]*u[j3];
              int_gauss += term;
            }
          }
        }
      }
    }
  }
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
  ofile.close();

  return 0;
}  // end of main program

// Help function for func_6D.
// It calculates the reciprocal value of the distance between r_1 and r_2.
double abs_dist(double x1, double y1, double z1, double x2, double y2, double z2) {
  double x1_x2 = (x1-x2)*(x1-x2);
  double y1_y2 = (y1-y2)*(y1-y2);
  double z1_z2 = (z1-z2)*(z1-z2);
  double r1_r2 = sqrt( x1_x2 + y1_y2 + z1_z2 );
  return r1_r2;
}

// Six dimensional integrand to be integrated.
// We simply exclude the integration points where the integrand diverges.
double func_6D(double x1, double y1, double z1, double x2, double y2, double z2) {
  double alpha = 2.0;
  double f_val = 0.0;
  double r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
  double r2 = sqrt(x2*x2 + y2*y2 + z2*z2);

  double r1_r2 = abs_dist(x1,y1,z1,x2,y2,z2);
  // Simply skip the parts where the integrand is singular:
  if (r1_r2 >= epsilon) { f_val = exp(-2.0*alpha*(r1 + r2)) / r1_r2; }
  return f_val;
}

/*
** The function
**              gauleg()
** takes the lower and upper limits of integration x1, x2, calculates
** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
** of length n of the Gauss--Legendre n--point quadrature formulae.
*/
void gauleg(double x1, double x2, double x[], double w[], int n)
{
  int         m,j,i;
  double      z1,z,xm,xl,pp,p3,p2,p1;
  double      const  pi = 3.14159265359;
  double      *x_low, *x_high, *w_low, *w_high;
  m  = (n + 1)/2;                             // roots are symmetric in the interval
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  x_low  = x;                                       // pointer initialization
  x_high = x + n - 1;
  w_low  = w;
  w_high = w + n - 1;

  for(i = 1; i <= m; i++) {                             // loops over desired roots
    z = cos(pi * (i - 0.25)/(n + 0.5));
      /*
      ** Starting with the above approximation to the ith root
      ** we enter the main loop of refinement by Newtons method.
      */
      do {
        p1 =1.0;
        p2 =0.0;
        /*
        ** loop up recurrence relation to get the
        ** Legendre polynomial evaluated at x
        */
        for(j = 1; j <= n; j++) {
          p3 = p2;
          p2 = p1;
          p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
        }
        /*
        ** p1 is now the desired Legrendre polynomial. Next compute
        ** ppp its derivative by standard relation involving also p2,
        ** polynomial of one lower order.
        */
        pp = n * (z * p1 - p2)/(z * z - 1.0);
        z1 = z;
        z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);
      /*
      ** Scale the root to the desired interval and put in its symmetric
      ** counterpart. Compute the weight and its symmetric counterpart
      */
      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
  }
} // End_ function gauleg()
