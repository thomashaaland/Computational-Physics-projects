#include <iostream>
#include "jacobi.h"

// This file runs a simple test on the jacobi method.
int main() {

  // Initialise the size of our matrices
  const int n = 2;

  // Initialise the A matrix
  double **A = new double*[n];
  for (int i = 0; i < n; i++) {
    A[i] = new double[n];
  }

  // Initialise the R matrix
  double **R = new double*[n];
  for (int i = 0; i < n; i++) {
    R[i] = new double[n];
  }

  // Fill A and R with values
  A[0][0] = 1; A[0][1] = 0.5; A[1][0] = -0.5; A[1][1] = 1;
  R[0][0] = 2; R[0][1] = 3; R[1][0] = -3; R[1][1] = 2;

  // Print A and R matrices
  std::cout << "A: " << std::endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "R: " << std::endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << R[i][j] << " ";
    }
    std::cout << "\n";
  }

  // Perform the jacobi method
  jacobi::jacobi_method(A, R, n);

  // Print the matrices again
  std::cout << "A: " << std::endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << A[i][j] << " ";
    }
    std::cout << "\n";
  }

  std::cout << "R: " << std::endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      std::cout << R[i][j] << " ";
    }
    std::cout << "\n";
  }
  
  for (int i = 0; i < n; i++) {
    delete[] A[i];
    delete[] R[i];
  }
  delete[] A;
  delete[] R;
  
  return 0;
}
