#include "jacobi.h"

namespace jacobi {
  void jacobi_method ( double ** A, double ** R, int n ) {
    //setting up the eigenvector matrix
    for ( int i = 0 ; i < n; i++) {
      for ( int j = 0 ; j < n; j++ ) {
	if ( i == j ) {
	  R[i][j] = 1.0;
	} else {
	  R[i][j] = 0.0;
	}
      }
    }
    
    int k, l;
    double epsilon = 1.0e-8;
    double max_number_terations = (double) n * (double) n * (double) n;
    int iterations = 0;
    double max_offdiag = jacobi::maxoffdiag ( A, &k, &l, n);

    while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_terations ) {
      max_offdiag = jacobi::maxoffdiag (A, &k, &l, n);
      jacobi::rotate (A, R, k, l, n);
      iterations++;
    }
    std::cout << "Number of iterations: " << iterations << "\n";
    return;
  }
  
  // Finding the maximum matrix element:
  double maxoffdiag ( double ** A, int * k, int * l, int n ) {
    double max = 0.0;
    
    for ( int i = 0; i < n; i++) {
      for ( int j = i +1; j < n; j++) {
	if ( fabs(A[i][j]) > max ) {
	  max = fabs(A[i][j]);
	  *l = i;
	  *k = j;
	}
      }
    }
    return max;
  }
  
  // Finding the values of cos and sin
  void rotate ( double ** A, double ** R, int k, int l, int n ) {
    double s, c;
    if ( A[k][l] != 0.0 ) {
      double t, tau;
      tau = (A[l][l] - A[k][k])/(2*A[k][l]);
      if ( tau > 0 ) {
	t = 1.0/(tau + sqrt(1.0 + tau*tau));
      } else {
	t = -1.0/( -tau + sqrt(1.0 + tau*tau));
      }
      
      c = 1/sqrt(1+t*t);
      s = c*t;
    } else {
      c = 1.0;
      s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    // changing the matrix elements with indices k and l
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0;
    A[l][k] = 0.0;
    // change elements
    for ( int i = 0; i < n; i++) {
      if ( i != k && i != l ) {
	a_ik = A[i][k];
	a_il = A[i][l];
	A[i][k] = c*a_ik - s*a_il;
	A[k][i] = A[i][k];
	A[i][l] = c*a_il + s*a_ik;
	A[l][i] = A[i][l];
      }
      // eigenvectors
      r_ik = R[i][k];
      r_il = R[i][l];
      R[i][k] = c*r_ik - s*r_il;
      R[i][l] = c*r_il + s*r_ik;
    }
    return;
  }
}
