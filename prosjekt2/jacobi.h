#pragma once
#include <armadillo>
#include <iostream>
#include <cmath>
#include "test.h"
namespace jacobi {
  // declare functions
  double Mmax (arma::mat &B, int &k,int &l, int n);
  void rotate(arma::mat &B, int k, int l, int n);
  arma::vec jacobi_eigenvalues(arma::mat &B, double epsilon, int maxit, int &k, int &l, int n, int &iterations);

}
