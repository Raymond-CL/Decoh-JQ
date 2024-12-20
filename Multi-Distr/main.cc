#include <gsl/gsl_integration.h>
#include <math.h>

#include <iostream>
#include <fstream>

const double Nc = 3;
const double Nf = 3;
const double b = 11.0 / 3.0 * Nc - 2.0 / 3.0 * Nf;
const double Lqcd = 0.25;
const double Q0 = 0.5;
const double lambda = log(Q0 / Lqcd);
const double ymin = 0.0;
const double ymax = 10.0;
const int Nn = 2;
const int Ny = 21;
double Pny_table[Nn+1][Ny];
double Sny_table[Nn+1][Ny];

struct input_params {
  int n_now;
  double y_now;
};

void initialize();
void print_Pny(std::ostream &o_obj);
inline double gamma(double yp) { return 4.0 * Nc / b / (yp + lambda); }

double get_Pny(int n, double y) {
  int i = 0;
  while (i < Ny - 1 && y > Pny_table[0][i]) {
    i++;
  }
  if (i == 0) {
    return Pny_table[n][0];
  } else if (i == Ny - 1) {
    return Pny_table[n][Ny - 1];
  } else {
    double x1 = Pny_table[0][i - 1];
    double x2 = Pny_table[0][i];
    double y1 = Pny_table[n][i - 1];
    double y2 = Pny_table[n][i];
    return y1 + (y - x1) * (y2 - y1) / (x2 - x1);
  }
}

double Sny_int(double yp, void *params) {
  struct input_params *p = (struct input_params *)params;
  int n = p->n_now;
  double y = p->y_now;
  if (n == 1) {
    return (y - yp) * gamma(yp);
  } else {
    return (y - yp) * gamma(yp) * get_Pny(n - 1, y);
  }
}

void set_Sny(int n) {
  // prepare integrtion workspace
  int npt = 1000;
  double result, error;
  double abserr = 1.0e-8;
  double relerr = 1.0e-8;
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(npt);
  gsl_function F;
  F.function = &Sny_int;

  for (int i = 0; i < Ny; i++) {
    double ypmin = ymin;
    double ypmax = Pny_table[n - 1][i];
    struct input_params ip = {n, ypmax};
    F.params = &ip;
    gsl_integration_qag(&F, ypmin, ypmax, abserr, relerr, npt,
                        GSL_INTEG_GAUSS15, w, &result, &error);
    Sny_table[n - 1][i] = result;
  }

  gsl_integration_workspace_free(w);
}

void set_Pny(int n) {
  if (n == 1) {
    for (int i = 0; i < Ny; i++) {
      Pny_table[n][i] = exp(-Sny_table[n - 1][i]);
    }
  } else {
    for (int i = 0; i < Ny; i++) {
      for (int j = 1; j < n; j++) {
        Pny_table[n][i] += j/(n-1)*Pny_table[n-j][i]*Sny_table[j][i];
      }
    }
  }
}


int main(void) {
  // initialize y values
  initialize();

  // run table
  for(int i=1;i<=Nn;i++){
    set_Sny(i);
    set_Pny(i);
  }

  // print table
  std::ofstream out_file("Pny_table.dat", std::ios::out);
  // print_Pny(std::cout);
  print_Pny(out_file);
  out_file.close();

  return 0;
}

void initialize() {
  for (int i = 0; i < Ny; i++) {
    Pny_table[0][i] = i * ymax / (Ny - 1);
  }
}

void print_Pny(std::ostream &o_obj) {
  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j <= Nn; j++) o_obj << Pny_table[j][i] << " \t ";
    o_obj << std::endl;
  }
}