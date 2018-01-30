// Copyright 2016 Alessandro Fabbri, Stefano Sinigardi

/***************************************************************************
This file is part of Inertial Analysis.

Inertial Analysis is free software : you can redistribute it and / or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Inertial Analysis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Inertial Analysis. If not, see <http://www.gnu.org/licenses/>.
***************************************************************************/

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "math_func.h"
#include "params.h"

//////// ANALYTICAL PARAMETRIZATION
double T0, T1, a, R, v1, x1, yy1, omega, period, T2, x2, y2, T3, x3, y3, T4, dt;
int N_SAMPLE;

void init_params()
{
  T0 = 0.0; // s
  T1 = T0 + 1.0; // s

  a = 1.0; // m/s^2
  R = 2.0; // m
  v1 = a * (T1 - T0); // m/s
  x1 = 0.5 * a * (T1 - T0) * (T1 - T0); // m
  yy1 = R; // m
  omega = v1 / R; // rad/s
  period = 2 * M_PI / omega; // s

  T2 = T1 + period / 8.; // s

  x2 = x1 + sqrt(2) * R; // m
  y2 = -R * (sqrt(2) - 1); // m

  T3 = T2 + period / 8.; // s

  x3 = x2; // m
  y3 = y2 + R; // m

  T4 = T3 + 1.0; // s

  N_SAMPLE = 2000;
  dt = (T4 - T0) / (N_SAMPLE - 1.); // s
}

double x(double t)
{
  double x;
  if (t <= T0) {
    x = 0.;
  } else if (t <= T1) {
    x = 0.5 * a * (t - T0) * (t - T0);
  } else if (t <= T2) {
    x = x1 + R * sin(omega * (t - T1));
  } else if (t <= T3) {
    x = x2 + R * cos(omega * (t - T2) - 3 / 4. * M_PI);
  } else if (t <= T4) {
    x = x3 + v1 * (t - T3) - 0.5 * a * (t - T3) * (t - T3);
  } else {
    x = x3 + 0.5;
  }

  return x;
}

double y(double t)
{
  double y;
  if (t <= T0) {
    y = 0.;
  } else if (t <= T1) {
    y = 0.;
  } else if (t <= T2) {
    y = yy1 - R * cos(omega * (t - T1));
  } else if (t <= T3) {
    y = y2 - R * sin(omega * (t - T2) - 3 / 4. * M_PI);
  } else if (t <= T4) {
    y = y3;
  } else {
    y = y3;
  }

  return y;
}

//////// MAIN
#define MAJOR 2
#define MINOR 0

int main(int argc, char** argv)
{
  std::cout << "Data Generator v" << MAJOR << "." << MINOR << std::endl;

  std::vector<double> t, a_x, a_y, v_x, v_y, r_x, r_y;
  std::vector<double> theta, Omega;
  std::complex<double> IU(0, 1); // Imaginary Unit

  init_params();

  std::cout << "Times (s) : " << T0 << " - " << T1 << " - " << T2 << " - " << T3 << " - " << T4 << std::endl;
  std::cout << "N_sample : " << N_SAMPLE << "\tdt (s) : " << dt << std::endl;

  for (int i = 0; i < N_SAMPLE; i++) {
    t.push_back(T0 + i * dt);
    r_x.push_back(x(t[i]));
    r_y.push_back(y(t[i]));
    v_x.push_back(0.0);
    v_y.push_back(0.0);
    a_x.push_back(0.0);
    a_y.push_back(0.0);
    theta.push_back(0.0);
    Omega.push_back(0.0);
  }

  v_x = forward_derivative(r_x, t);
  v_y = forward_derivative(r_y, t);

  a_x = forward_derivative(v_x, t);
  a_y = forward_derivative(v_y, t);

  for (size_t i = 0; i < t.size(); i++) {
    double th = atan2(v_y[i], v_x[i]);
    if (th < 0)
      th += 2 * M_PI;
    theta[i] = th;
  }
  Omega = forward_derivative(theta, t);

  char out_name[20];
  sprintf(out_name, "data_lab.txt");
  FILE* out_data = fopen(out_name, "w");
  fprintf(out_data, "#%5s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "Index", "t", "r_x", "r_y", "v_x", "v_y", "a_x", "a_y", "theta", "omega");
  for (size_t i = 0; i < t.size(); i++)
    fprintf(out_data, "%6zu %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n",
        i, t[i], r_x[i], r_y[i], v_x[i], v_y[i], a_x[i], a_y[i], theta[i], Omega[i]);
  fclose(out_data);

  std::vector<std::complex<double>> a_loc, a_fix;
  std::vector<double> a_t, a_n;
  for (size_t i = 0; i < t.size(); i++) {
    a_fix.push_back(std::complex<double>(a_x[i], a_y[i]));
    a_loc.push_back(exp(-1. * IU * theta[i]) * a_fix[i]);
    a_t.push_back(a_loc[i].real());
    a_n.push_back(a_loc[i].imag());
  }

  out_data = fopen("data_loc.txt", "w");
  fprintf(out_data, "#%5s %16s %16s %16s %16s \n", "i", "t", "a_t", "a_n", "omega_z");
  for (size_t i = 0; i < t.size(); i++)
    fprintf(out_data, "%6zu %16.10lf %16.10lf %16.10lf %16.10lf\n",
        i, t[i], a_t[i], a_n[i], Omega[i]);
  fclose(out_data);

  std::default_random_engine generator;
  std::normal_distribution<double> noise_acc(0.0, 2.5);
  std::normal_distribution<double> noise_gyr(0.0, 0.85);

  out_data = fopen("data_inertial.txt", "w");
  fprintf(out_data, "#%15s %16s %16s %16s %16s %16s %16s %16s \n", "t", "|v|", "a_t", "a_n", "a_z", "omega_x", "omega_y", "omega_z");
  for (size_t i = 0; i < t.size(); i++)
    fprintf(out_data, "%16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf %16.10lf\n",
        t[i], 0.0, a_t[i] + noise_acc(generator), a_n[i] + noise_acc(generator), 0.0, 0.0, 0.0, Omega[i] + noise_gyr(generator));
  fclose(out_data);

  return 0;
}
