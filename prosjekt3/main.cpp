#include "gauleg.h"
//   Project 3 a)
//   A first approach to calculating the 6D integral numerically
//   with Gauss-Legendre quadrature.

//using namespace std;
std::ofstream ofile;

//   Main function begins here
int main()
{
  // Initial read in of some numbers:
  int n;
  int k = 50;
  int A[k];
  for (int i = 1; i <= k; i++) {A[i-1] = i;}

  // Values of n to be tested 
  int size_of_A = sizeof(A)/sizeof(A[0]);
  
  double x1, x2; // Limits of the integral in Cartesian coordinates
  double const  pi = 3.14159265359;
  double exact_integral = 5*pi*pi/256; // Exact value of the integral to be evaluated
  double int_gauss = 0.0;
  double relatative_error = 0.0;
  double calculation_time = 0.0;

  //gauleg::func_inspect(1);
  
  /*
  // Inspect the chosen function
  std::string inspect_path = "func_inspect.csv";
  ofile.open(inspect_path);
  //gauleg::func_inspect(gauleg::func_6D);
  ofile << "X, val" << std::endl;
  double X;
  std::vector<double> args;
  args.reserve(6);
  for (int i = 0; i < 101; i++) {
    X = 2*((double) i) / 101 - 1;
    //args = {X, 0, 0, 0, 0, 0};
    ofile << X << ", " << gauleg::func_6D(X, 0, 0, 0, 0, 0) << std::endl;
  }
  ofile.close();
  */
  
  std::string outfilename = "Gauss_Legendre_table.txt";
  //   reserve space in memory for vectors containing the mesh points
  //   weights and function values for the use of the gauss-legendre method
  double *x;
  double *w;

  std::cout << "Read in integration limits in the format: lower upper" << std::endl;
  std::cin >> x1 >> x2;
  std::cout << "Got the limits " << x1 << " " << x2 << std::endl;

  std::cout << "Creating file " << outfilename << std::endl;
  ofile.open(outfilename);
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << " Integration limits in each dimension: x1 = " << x1 << " and x2 = " << x2 << std::endl;
  ofile << " n,        Result with Gauss-Legendre,      Exact result,       Relative error,     Calculation time [s]" << std::endl;
  std::cout << "File created. Looping over " << size_of_A << " values" << std::endl;
  // loop over all values of n:
  for (int l=0; l < size_of_A; l++) {
    //inits
    std::cout << "\r" << (int)((double)(l+1)/(size_of_A)*100) << "% complete" << std::flush;
    n = A[l];
    x = new double [n];
    w = new double [n];
    int_gauss = 0.;
    calculation_time = 0.0;
    clock_t start, finish;
    start = clock();

    
    // trying the 6D integral with Gauss-Legendre method:
    int_gauss = gauleg::gauss_legendre(x1, x2, x, w, n);

    finish = clock();
    calculation_time = (finish - start)/(double)CLOCKS_PER_SEC;
    relatative_error = fabs(int_gauss - exact_integral)/exact_integral;
    
    ofile << std::setw(5) << std::setprecision(5) << n << ",";
    ofile << std::setw(25) << std::setprecision(10) << int_gauss << ",";
    ofile << std::setw(25) << std::setprecision(10) << exact_integral << ",";
    ofile << std::setw(20) << std::setprecision(4) << relatative_error << ",";
    ofile << std::setw(20) << std::setprecision(4) << calculation_time << std::endl;
    delete [] x;
    delete [] w;
  }
  std::cout << std::endl;

  ofile.close();
  
  return 0;
}  // end of main program


