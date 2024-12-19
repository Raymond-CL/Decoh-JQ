#include <gsl/gsl_integration.h>
#include <math.h>

#include <iostream>

const double Nc = 3;
const double Nf = 4;
const double b = 11.0 / 3.0 * Nc - 2.0 / 3.0 * Nf;
const double Lqcd = 0.25;
const double Q0 = 0.5;
const double lambda = log(Q0 / Lqcd);
const double y = 10.0;

const int nmax = 3;
const int ny = 11;
double Pn_table[nmax][ny];

inline double gamma(double yp) { return 4.0 * Nc / b / (yp + lambda); }

double f(double yp, void* params) { return (y - yp) * gamma(yp); }

double P1(double y) { return 0; }

int main(void) {
  for (int i = 0; i < ny; i++) {
    Pn_table[0][i] = i * y / (ny - 1);
    // std::cout << Pn_table[0][i] << std::endl;
  }

  for (int i = 0; i < ny; i++) {
    Pn_table[1][i] = P1(Pn_table[0][i]);
  }

  for (int i = 0; i < ny; i++) {
    std::cout << Pn_table[0][i] << '\t' << Pn_table[1][i] << std::endl;
  }

  size_t n = 100;
  double ypmin = 0.0;
  double ypmax = y;
  double epsabs = 1.0;
  double epsrel = 1.0;
  double result, error;
  double expected =
      -4.0 * Nc * (y + (y + lambda) * log(lambda / (y + lambda))) / b;

  gsl_integration_workspace* w = gsl_integration_workspace_alloc(n);

  gsl_function F;
  F.function = &f;

  gsl_integration_qag(&F, ypmin, ypmax, epsabs, epsrel, n, GSL_INTEG_GAUSS21, w,
                      &result, &error);

  printf("result          = % .18f\n", result);
  printf("exact result    = % .18f\n", expected);
  printf("estimated error = % .18f\n", error);
  printf("actual error    = % .18f\n", result - expected);
  printf("intervals       = %zu\n", w->size);

  gsl_integration_workspace_free(w);

  return 0;
}