#pragma once

// Declare functions:
namespace jacobi {
  void jacobi_method(double ** A, double ** R, int n);
  double maxoffdiag(double ** A, int * k, int * l, int n);
  void rotate ( double ** A, double ** R, int k, int l, int n );
}
