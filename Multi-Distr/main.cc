#include <gsl/gsl_integration.h>
#include <math.h>

#include <iostream>

const double ymax = 10.0;
const int nmax = 5;
const int ny = 51;
double Pn_table[nmax][ny];
double Sn_table[nmax][ny];

inline double gamma(double yp) {
  const double Nc = 3;
  const double Nf = 3;
  const double b = 11.0 / 3.0 * Nc - 2.0 / 3.0 * Nf;
  const double Lqcd = 0.25;
  const double Q0 = 0.5;
  const double lambda = log(Q0 / Lqcd);
  return 4.0 * Nc / b / (yp + lambda);
}

double f(double yp, void* params) {
  double y = *(double*)params;
  return (y - yp) * gamma(yp);
}

double P1(double y) {
  size_t n = 10;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(n);
  double ypmin = 0.0;
  double ypmax = y;
  double epsabs = 0.0;
  double epsrel = 0.1;
  double result, error;
  gsl_function F;
  F.function = &f;
  F.params = &y;
  gsl_integration_qag(&F, ypmin, ypmax, epsabs, epsrel, n, GSL_INTEG_GAUSS15, w,
                      &result, &error);
  gsl_integration_workspace_free(w);
  return exp(-result);
}

double linear_interpolate(const double* x, const double* y, int size,
                          double x_value) {
  if (size < 2) {
    throw std::invalid_argument(
        "Array size must be at least 2 for interpolation.");
  }

  // Handle the case where x_value is out of bounds
  if (x_value <= x[0]) {
    return y[0];
  }
  if (x_value >= x[size - 1]) {
    return y[size - 1];
  }

  // Find the interval [x[i], x[i+1]] that contains x_value
  int i = 0;
  while (x_value > x[i + 1] && i < size - 2) {
    i++;
  }

  // Perform linear interpolation
  double x0 = x[i];
  double x1 = x[i + 1];
  double y0 = y[i];
  double y1 = y[i + 1];

  return y0 + (y1 - y0) * (x_value - x0) / (x1 - x0);
}

double g(double yp, void* params) {
  double y = *(double*)params;
  return (y - yp) * gamma(yp) *
         linear_interpolate(Pn_table[0], Pn_table[1], ny, yp);
}

double S1(double y) {
  size_t n = 10;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(n);
  double ypmin = 0.0;
  double ypmax = y;
  double epsabs = 0.0;
  double epsrel = 0.1;
  double result, error;
  gsl_function F;
  F.function = &g;
  F.params = &y;
  gsl_integration_qag(&F, ypmin, ypmax, epsabs, epsrel, n, GSL_INTEG_GAUSS15, w,
                      &result, &error);
  gsl_integration_workspace_free(w);
  return result;
}

int main(void) {
  for (int i = 0; i < ny; i++) {
    Pn_table[0][i] = i * ymax / (ny - 1);
  }

  for (int i = 0; i < ny; i++) {
    Pn_table[1][i] = P1(Pn_table[0][i]);
  }

  for (int i = 0; i < ny; i++) {
    Pn_table[2][i] = S1(Pn_table[0][i]);
  }

  for (int i = 0; i < ny; i++) {
    Pn_table[3][i] = Pn_table[1][i] * Pn_table[2][i];
  }

  for (int i = 0; i < ny; i++) {
    std::cout << Pn_table[0][i] << '\t' << Pn_table[1][i] << '\t'
              << Pn_table[3][i] << std::endl;
  }

  return 0;
}