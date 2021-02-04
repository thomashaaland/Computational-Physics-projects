#include "jacobi.h"

namespace jacobi {
  arma::vec jacobi_eigenvalues(arma::mat &B, double epsilon, int maxit, int &k, int &l, int n, int &iterations) {
    // Rotate the Bmatrix to find eigenvalues
    iterations = 0;
    double Bmax = 10;
    while (abs(Bmax) > epsilon && iterations < maxit) {
      // find k, l max value
      Bmax = Mmax(B,k,l,n);
      rotate(B, k, l, n);
      iterations++;
      
    }
    arma::vec eigenB = sort(B.diag());
    return eigenB;
  }
  
  void rotate(arma::mat &B, int k, int l, int n) {
    double a_kk, a_ll, a_ik, a_il;
    double c, s, t, tau;
    
    // First, find the value tau = cot(2x), and for t
    if (B(k,l) != 0.0) {
      tau = (B(l,l)-B(k,k))/(2*B(k,l));
      if (tau > 0.0) {
	t = 1.0/(tau+sqrt(1+tau*tau));
      } else {
	t = 1.0/(tau-sqrt(1+tau*tau));
      }
      c = 1/sqrt(1+t*t);
      s = t*c;
    } else {
      c = 1.0;
      s = 0.0;
    }
    a_kk = B(k,k);
    a_ll = B(l,l);
    // changing the matrix elements with indices k and l
    B(k,k) = c*c*a_kk - 2.0*c*s*B(k,l) + s*s*a_ll;
    B(l,l) = s*s*a_kk + 2.0*c*s*B(k,l) + c*c*a_ll;
    B(k,l) = 0.0; // hard-coding of the arma::zeros
    B(l,k) = 0.0;
    
    // and then we change the remaining elements
    for ( int i = 0; i < n; i++ ) {
      if ( i != k && i != l ) {
	a_ik = B(i,k);
	a_il = B(i,l);
	B(i,k) = c*a_ik - s*a_il;
	B(k,i) = B(i,k);
	B(i,l) = c*a_il + s*a_ik;
	B(l,i) = B(i,l);
      }
    }
  }
  
  // Find the largest offdiagonal value in an symmetric matrix and return k and l coordinates
  
  double Mmax (arma::mat &B, int &k,int &l, int n){
    int i, j;
    double Bmax=0.;
    for(i = 0 ; i < n ; i++){
      for(j = i+1 ; j < n ; j++) {
	if (abs(Bmax) < abs(B(i,j))) {
	  Bmax = B(i,j);
	  k = i;
	  l = j;
	}
      }
    }
    return Bmax;
  }
  
}
