#pragma once
#include <thread>
#include <cmath>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


// Definitions
#define EPS 3.0e-14
#define MAXIT 10
#define ZERO 1.0E-10
#define epsilon 1.0E-6

namespace gauleg {
  //     Here we define various functions called by the main program:
  void gauss_laguerre(double *, double *, int, double);
  void gauleg(double, double, double *, double *, int);
  double func_6D(double, double, double, double, double, double);
  double abs_dist(double, double, double, double, double, double);
  double gauss_legendre(double a, double b, double *r, double *u, int n);
}
