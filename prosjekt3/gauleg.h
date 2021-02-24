#pragma once
#include <cmath>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <time.h>
#include <vector>

// Definitions
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define epsilon 1.0E-6

namespace gauleg {
  //     Here we define various functions called by the main program:
  void gauleg(double, double, double *, double *, int);
  double func_6D(double, double, double, double, double, double);
  double abs_dist(double, double, double, double, double, double);
  double gauss_legendre(double x1, double x2, double *x, double *w, int n);
  
}
