#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gauleg.h"

// Help function for func_6D.
// It calculates the reciprocal value of the distance between r_1 and r_2.
double gauleg::abs_dist(double x1, double y1, double z1, double x2, double y2, double z2) {
  double x1_x2 = (x1-x2)*(x1-x2);
  double y1_y2 = (y1-y2)*(y1-y2);
  double z1_z2 = (z1-z2)*(z1-z2);
  double r1_r2 = sqrt( x1_x2 + y1_y2 + z1_z2 );
  return r1_r2;
}

// Six dimensional integrand to be integrated.
// We simply exclude the integration points where the integrand diverges.
double gauleg::func_6D(double x1, double y1, double z1, double x2, double y2, double z2) {
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
void gauleg::gauleg(double x1, double x2, double x[], double w[], int n)
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
