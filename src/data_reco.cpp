// Copyright 2016 Alessandro Fabbri, Stefano Sinigardi

/***************************************************************************
This file is part of inertial_tools.

inertial_tools is free software : you can redistribute it and / or modify
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
#include <fstream>
#include <iostream>
#include <vector>

// boost include
#include <boost/algorithm/string.hpp>
#include <boost/utility.hpp>

// jsoncons include
#include "jsoncons/json.hpp"

#include "io_lib.hpp"
#include "math_func.h"
#include "params.h"

#define MAJOR_VERSION 2
#define MINOR_VERSION 1

void usage(char* progname)
{
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " init_file data_file" << std::endl;
  std::cout << "       init_file must be a valid initial condition JSON file (launch with no arguments to create a template)" << std::endl;
  std::cout << "       data_file must be \"inertial\" PHYSYCOM standard compliant" << std::endl
            << std::endl;
  std::cout << "Usage: " << progname << " -init_t" << std::endl;
  std::cout << "       create a template init JSON file" << std::endl;
  exit(-3);
}

int main(int argc, char** argv)
{
  std::cout << "Trajectory Reconstruction 2D v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl
            << std::endl;

  std::string init_file, data_file;
  if (argc == 3) {
    init_file = argv[1];
    data_file = argv[2];
  } else if (argc == 2 && std::string(argv[1]) == "-init_t") {
    std::cout << "Generating initial condition JSON template \"init.json\"" << std::endl
              << std::endl;
    std::ofstream init_json("init.json");
    init_json << R"({
  "r0_x" : 0.0,
  "r0_y" : 0.0,
  "theta0" : 0.0,
  "v0_x" : 0.0,
  "v0_y" : 0.0,
  "dt"   : 1e-3, // [optional] use it to force the timestamp to be exactly equally spaced
  "acc_conversion" : [ 9.81, 9.81, 9.81],
  "gyr_conversion" : [ 0.017453, 0.017453, 0.017453],
  "speed_conversion" : 0.2778
}
)";
    init_json.close();
    exit(-1);
  } else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
  }

  std::cout << "Reconstructing trajectory for input : " << data_file << std::endl;

  std::vector<std::complex<double>> a_fixed, a_local;
  std::vector<double> timestep, omega_z, a_local_x, a_local_y;
  std::vector<double> theta, a_x, a_y, v_x, v_y, r_x, r_y;
  double theta0, v0_x, v0_y, r0_x, r0_y, dt = 0.0;
  double acc_conversion[3], gyr_conversion[3], speed_conversion;

  // Imaginary Unit
  std::complex<double> IU(0, 1);

  // parsing initial condition file
  jsoncons::json init;
  try {
    init = jsoncons::json::parse_file(init_file);
  } catch (std::exception& e) {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    exit(-7);
  }
  r0_x = init.has_member("r0_x") ? init["r0_x"].as<double>() : 0.0;
  r0_y = init.has_member("r0_y") ? init["r0_y"].as<double>() : 0.0;
  theta0 = init.has_member("theta0") ? init["theta0"].as<double>() : 0.0;
  v0_x = init.has_member("v0_x") ? init["v0_x"].as<double>() : 0.0;
  v0_y = init.has_member("v0_y") ? init["v0_y"].as<double>() : 0.0;
  dt = init.has_member("dt") ? init["dt"].as<double>() : 1.e-3;
  speed_conversion = init.has_member("speed_conversion") ? init["speed_conversion"].as<double>() : 1 / 3.6;
  for (int i = 0; i < 3; i++) {
    acc_conversion[i] = init.has_member("acc_conversion") ? init["acc_conversion"][i].as<double>() : GRAV;
    gyr_conversion[i] = init.has_member("gyr_conversion") ? init["gyr_conversion"][i].as<double>() : 1 / RAD_TO_DEG;
  }

  std::cout << "Inital condition : " << data_file << std::endl
            << "r0 = ( " << r0_x << " , " << r0_y << " )" << std::endl
            << "v0 = ( " << v0_x << " , " << v0_y << " )" << std::endl
            << "theta0 = " << theta0 << std::endl
            << "dt = " << (init.has_member("dt") ? std::to_string(dt) : "dynamic") << std::endl;
  std::cout << "Conversion factors : " << std::endl;
  std::cout << "Speed " << speed_conversion << std::endl;
  std::cout << "Acc [ ";
  for (int i = 0; i < 3; i++)
    std::cout << acc_conversion[i] << ((i == 2) ? " ]" : " , ");
  std::cout << std::endl;
  std::cout << "Gyr [ ";
  for (int i = 0; i < 3; i++)
    std::cout << gyr_conversion[i] << ((i == 2) ? " ]" : " , ");
  std::cout << std::endl;

  // parsing data file
  std::vector<std::vector<std::string>> file_tokens;
  file_tokens = Read_from_file(data_file);
  double timenow = -dt;
  for (size_t i = 0; i < file_tokens.size(); i++) {
    if (init.has_member("dt")) {
      timenow += dt;
      timestep.push_back(timenow);
    } else {
      timestep.push_back(atof(file_tokens[i][TIMESTAMP_INDEX].c_str()));
    }
    omega_z.push_back(atof(file_tokens[i][GZ_INDEX].c_str()));
    a_local_x.push_back(atof(file_tokens[i][AX_INDEX].c_str()));
    a_local_y.push_back(atof(file_tokens[i][AY_INDEX].c_str()));
  }

  //////// PHYSICS
  // converting data to SI units
  std::transform(omega_z.begin(), omega_z.end(), omega_z.begin(), bind1st(std::multiplies<double>(), gyr_conversion[2])); // omega_z MUST be in RADIANS
  std::transform(a_local_x.begin(), a_local_x.end(), a_local_x.begin(), bind1st(std::multiplies<double>(), acc_conversion[0])); // a_x MUST be in m/s^2
  std::transform(a_local_y.begin(), a_local_y.end(), a_local_y.begin(), bind1st(std::multiplies<double>(), acc_conversion[1])); // a_y MUST be in m/s^2

  // integrators
  theta = Integrate(timestep, omega_z, theta0);

  for (size_t i = 0; i < a_local_x.size(); i++) {
    a_local.push_back(std::complex<double>(a_local_x[i], a_local_y[i]));
    a_fixed.push_back(exp(1. * IU * theta[i]) * a_local[i]);
    a_x.push_back(a_fixed[i].real());
    a_y.push_back(a_fixed[i].imag());
  }

  v_x = Integrate(timestep, a_x, v0_x);
  v_y = Integrate(timestep, a_y, v0_y);

  r_x = Integrate(timestep, v_x, r0_x);
  r_y = Integrate(timestep, v_y, r0_y);

  //////// Saving data to file
  if (data_file.substr(0, 2) == ".\\")
    data_file = data_file.substr(2, data_file.size()); // remove leading ".\" in filename, if any
  std::string output = data_file.substr(0, data_file.size() - 4) + "_trajectory.txt";
  FILE* out_data = fopen(output.c_str(), "w");
  fprintf(out_data, "#%5s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n", "Index", "Time", "r_x", "r_y", "v_x", "v_y", "a_x", "a_y", "theta_z", "w_z");
  for (size_t i = 0; i < r_x.size(); i++)
    fprintf(out_data, "%6zu %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n",
        i, timestep[i] - timestep[0], r_x[i], r_y[i], v_x[i], v_y[i], a_x[i], a_y[i], theta[i], omega_z[i]);
  fclose(out_data);
  std::cout << "Created output data file : " << output << std::endl;

  //////// Generating gnuplot scripts
  std::string plot_basename = data_file.substr(0, data_file.size() - 4) + "_plot";
  std::ofstream gnuplot(plot_basename + ".plt");
  gnuplot << R"(#!/usr/bin/gnuplot
reset
set terminal pngcairo dashed size 900, 700 enhanced font 'Verdana,10'
set output ')"
          << plot_basename << R"(.png'
set size ratio 1
set multiplot layout 2, 2 title "Trajectory reconstruction" font ",14"
# Styles
linew = 1.2
points = 0.3
set style line  21 lc rgb '#0072bd' pointtype 7 pointsize points  # blue
set style line  22 lc rgb '#d95319' pointtype 7 pointsize points  # orange
set style line  23 lc rgb '#77ac30' pointtype 7 pointsize points  # green
set style line  24 lc rgb '#a2142f' pointtype 7 pointsize points  # red
set style line 101 lc rgb '#000000' lt 1 lw 1                     # black
set style line 102 lc rgb '#d6d7d9' lt 1 lw 1                     # gray
# Border xy
set border 3 front ls 101
set tics nomirror out scale 0.75
set format '%g'
set border linewidth 1.5
# Legend
set key top
# Plot
input = ")"
          << output << R"("
set size 0.5, 0.5
set origin -0.05, 0.45
set xlabel 't'
set ylabel 'x'
plot input using 2:3 notitle with points ls 21
set size 0.5, 0.5
set origin 0.35, 0.45
set xlabel 't'
set ylabel 'y'
plot input using 2:4 notitle with points ls 22
set size 0.5, 0.5
set origin -0.05, 0.0
set xlabel 'x'
set ylabel 'y'
plot input using 3:4 notitle with points ls 23
unset multiplot
)";
  std::cout << "Created gnuplot script file : " << plot_basename << ".plt" << std::endl;

  return 0;
}
