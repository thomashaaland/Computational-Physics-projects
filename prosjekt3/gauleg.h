#pragma once
#include <assert.h>
#include <cmath>
#include <cmath>
#include <fstream>
#include <functional>
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
  //double func_6D(std::vector<double> &args);
  double abs_dist(double, double, double, double, double, double);
  double gauss_legendre(double x1, double x2, double *x, double *w, int n);
  //template <typename ...ArgsT>
  //void func_inspect(double func(std::vector<double> &args));
  template <typename T, typename... ArgsT>
  void func_inspect(T arg1, ArgsT... arg2);
}
