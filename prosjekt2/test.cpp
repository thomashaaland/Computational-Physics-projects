#include "test.h"

namespace test {
  void UnitTests () {
    std::cout << "Test results: " << std::endl;
    int a = test::Mmax_test();
    int b = test::rotate_test();
    int c = test::jacobi_test();
    int tot = a+b+c;
    int number_of_tests = 3;
    std::cout << "Successful tests: " << tot << " out of "
	    << number_of_tests << " tests performed!" << std::endl;
  }
  
  int rotate_test(){
    int n = 2;
    arma::mat C = arma::ones(n,n);
    arma::vec eigen_test;
    jacobi::rotate(C,0,1,n);
    eigen_test = sort(C.diag());
    if (C(0,1) == 0.0 && C(1,0) == 0.0 && eigen_test(0) <= 1e-8 && eigen_test(1) - 2.0 <= 1e-8) {
      std::cout << "Rotate is working as advertised!" << std::endl;
      return 1;
    } else {
      std::cout << "Rotate is a failure!" << std::endl;
      return 0;
    }
  }
  
  int jacobi_test (){
    int n = 2;
    int l, k;
    int maxit = 100;
    int iterations;
    double epsilon = 1e-8;
    arma::mat C = arma::ones(n,n);
    arma::vec eigen_test;
    eigen_test = sort(jacobi::jacobi_eigenvalues(C, epsilon, maxit, k, l, n, iterations));
    if (eigen_test(0) <= epsilon && eigen_test(1) - 2.0 <= epsilon && iterations-1 == 1) {
      std::cout << "Jacobi_eigenvalues is working as advertised!" << std::endl;
      return 1;
    } else {
      std::cout << "Jacobi_eigenvalues is a failure!" << std::endl;
      return 0;
    }
  }
  
  int Mmax_test () {
    int n = 3;
    arma::mat C = arma::zeros<arma::mat>(n,n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	C(i,j) = n*i+j+1;
      }
    }
    int k, l;
    double Cmax = jacobi::Mmax(C, k, l, n);
    if (k == 1 && l == 2 && Cmax == 6){
      std::cout << "Mmax is working as advertised!" << std::endl;
      return 1;
    } else {
      std::cout << "Mmax is a failure!" << std::endl;
      return 0;
    }
  }
}
