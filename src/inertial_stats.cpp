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


#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif

#include "params.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <limits>

#include "io_lib.hpp"
#include "math_lib.h"


#define MAJOR_VERSION       0
#define MINOR_VERSION       2

void usage(char * progname) {
  std::cerr << "Usage: " << progname << " path/to/data/file" << std::endl;
}

int main(int argc, char **argv) {
  std::cout << "inertial_stats v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  // Command line parsing
  std::string input_file;
  if (argc > 1) {
    input_file = argv[1];
    if (input_file.substr(0, 2) == ".\\") input_file = input_file.substr(2, input_file.size() - 2);
    if (input_file.substr(0, 2) == "./") input_file = input_file.substr(2, input_file.size() - 2);
  }
  else {
    std::cerr << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
    exit(-1);
  }

  // Data parsing, convertion, storage
  vector< vector<string> > file_tokens = Read_from_file(input_file);
  vector< vector<double> > data = tokens_to_double(file_tokens);
  vector<double> ax, ay, az, gx, gy, gz;
  for (auto line : data) {
    ax.push_back(line[AX_INDEX]);
    ay.push_back(line[AY_INDEX]);
    az.push_back(line[AZ_INDEX]);
    gx.push_back(line[GX_INDEX]);
    gy.push_back(line[GY_INDEX]);
    gz.push_back(line[GZ_INDEX]);
  }

  // Ranges
  auto   ax_mm = std::minmax_element(ax.begin(), ax.end());
  double ax_min = *ax_mm.first;
  double ax_max = *ax_mm.second;

  auto   ay_mm = std::minmax_element(ay.begin(), ay.end());
  double ay_min = *ay_mm.first;
  double ay_max = *ay_mm.second;

  auto   az_mm = std::minmax_element(az.begin(), az.end());
  double az_min = *az_mm.first;
  double az_max = *az_mm.second;

  auto   gx_mm = std::minmax_element(gx.begin(), gx.end());
  double gx_min = *gx_mm.first;
  double gx_max = *gx_mm.second;

  auto   gy_mm = std::minmax_element(gy.begin(), gy.end());
  double gy_min = *gy_mm.first;
  double gy_max = *gy_mm.second;

  auto   gz_mm = std::minmax_element(gz.begin(), gz.end());
  double gz_min = *gz_mm.first;
  double gz_max = *gz_mm.second;

  // Evaluating statistics
  VEC3D ave_acc, ave_gyr;
  MAT3D quad_acc, quad_gyr, quad_mix;
  memset(&ave_acc, 0, sizeof(VEC3D));
  memset(&ave_gyr, 0, sizeof(VEC3D));
  memset(&quad_acc, 0, sizeof(MAT3D));
  memset(&quad_gyr, 0, sizeof(MAT3D));
  memset(&quad_mix, 0, sizeof(MAT3D));

  for (size_t i = 0; i < ax.size(); i++) {
    ave_acc.x += ax[i];
    ave_acc.y += ay[i];
    ave_acc.z += az[i];
    ave_gyr.x += gx[i];
    ave_gyr.y += gy[i];
    ave_gyr.z += gz[i];

    quad_acc.xx += ax[i] * ax[i]; quad_acc.xy += ax[i] * ay[i]; quad_acc.xz += ax[i] * az[i];
    quad_acc.yx += ay[i] * ax[i]; quad_acc.yy += ay[i] * ay[i]; quad_acc.xz += ay[i] * az[i];
    quad_acc.zx += az[i] * ax[i]; quad_acc.zy += az[i] * ay[i]; quad_acc.zz += az[i] * az[i];

    quad_gyr.xx += gx[i] * gx[i]; quad_gyr.xy += gx[i] * gy[i]; quad_gyr.xz += gx[i] * gz[i];
    quad_gyr.yx += gy[i] * gx[i]; quad_gyr.yy += gy[i] * gy[i]; quad_gyr.xz += gy[i] * gz[i];
    quad_gyr.zx += gz[i] * gx[i]; quad_gyr.zy += gz[i] * gy[i]; quad_gyr.zz += gz[i] * gz[i];

    quad_mix.xx += ax[i] * gx[i]; quad_mix.xy += ax[i] * gy[i]; quad_mix.xz += ax[i] * gz[i];
    quad_mix.yx += ay[i] * gx[i]; quad_mix.yy += ay[i] * gy[i]; quad_mix.xz += ay[i] * gz[i];
    quad_mix.zx += az[i] * gx[i]; quad_mix.zy += az[i] * gy[i]; quad_mix.zz += az[i] * gz[i];
  }
  multiply_vec3d(1. / math_float(ax.size()), &ave_acc);
  multiply_vec3d(1. / math_float(ax.size()), &ave_gyr);
  multiply_mat3d(1. / math_float(ax.size()), &quad_acc);
  multiply_mat3d(1. / math_float(ax.size()), &quad_gyr);
  multiply_mat3d(1. / math_float(ax.size()), &quad_mix);

  // Statistics estimator
  VEC3D devstd_acc, devstd_gyr;
  MAT3D cov_acc, cov_gyr, cov_mix;
  memset(&cov_acc, 0, sizeof(MAT3D));
  memset(&cov_gyr, 0, sizeof(MAT3D));
  memset(&cov_mix, 0, sizeof(MAT3D));

  devstd_acc.x = sqrt(quad_acc.xx - ave_acc.x*ave_acc.x);
  devstd_acc.y = sqrt(quad_acc.yy - ave_acc.y*ave_acc.y);
  devstd_acc.z = sqrt(quad_acc.zz - ave_acc.z*ave_acc.z);

  devstd_gyr.x = sqrt(quad_gyr.xx - ave_gyr.x*ave_gyr.x);
  devstd_gyr.y = sqrt(quad_gyr.yy - ave_gyr.y*ave_gyr.y);
  devstd_gyr.z = sqrt(quad_gyr.zz - ave_gyr.z*ave_gyr.z);

  cov_acc.xx = quad_acc.xx - ave_acc.x*ave_acc.x; cov_acc.xy = quad_acc.xy - ave_acc.x*ave_acc.y; cov_acc.xz = quad_acc.xz - ave_acc.x*ave_acc.z;
  cov_acc.xy = quad_acc.xy - ave_acc.x*ave_acc.y; cov_acc.yy = quad_acc.yy - ave_acc.y*ave_acc.y; cov_acc.yz = quad_acc.yz - ave_acc.y*ave_acc.z;
  cov_acc.xz = quad_acc.xz - ave_acc.x*ave_acc.z; cov_acc.zy = quad_acc.zy - ave_acc.z*ave_acc.y; cov_acc.zz = quad_acc.zz - ave_acc.z*ave_acc.z;

  cov_gyr.xx = quad_gyr.xx - ave_gyr.x*ave_gyr.x; cov_gyr.xy = quad_gyr.xy - ave_gyr.x*ave_gyr.y; cov_gyr.xz = quad_gyr.xz - ave_gyr.x*ave_gyr.z;
  cov_gyr.xy = quad_gyr.xy - ave_gyr.x*ave_gyr.y; cov_gyr.yy = quad_gyr.yy - ave_gyr.y*ave_gyr.y; cov_gyr.yz = quad_gyr.yz - ave_gyr.y*ave_gyr.z;
  cov_gyr.xz = quad_gyr.xz - ave_gyr.x*ave_gyr.z; cov_gyr.zy = quad_gyr.zy - ave_gyr.z*ave_gyr.y; cov_gyr.zz = quad_gyr.zz - ave_gyr.z*ave_gyr.z;

  // Output
  std::string out_filename = input_file.substr(0, input_file.size() - 4) + "_stats.txt";
  std::ofstream data_out(out_filename);
  data_out << "INERTIAL DATA size : " << ax.size() << std::endl;
  data_out << "ACC averages       : " << ave_acc.x << "  " << ave_acc.y << "  " << ave_acc.z << std::endl;
  data_out << "ACC std_dev        : " << devstd_acc.x << "  " << devstd_acc.y << "  " << devstd_acc.z << std::endl;
  data_out << "ACC estimators     : " << devstd_acc.x / fabs(ave_acc.x) << "  " << devstd_acc.y / fabs(ave_acc.y) << "  " << devstd_acc.z / fabs(ave_acc.z) << std::endl;
  data_out << "GYR averages       : " << ave_gyr.x << "  " << ave_gyr.y << "  " << ave_gyr.z << std::endl;
  data_out << "GYR std_dev        : " << devstd_gyr.x << "  " << devstd_gyr.y << "  " << devstd_gyr.z << std::endl;
  data_out << "GYR estimators     : " << devstd_gyr.x / fabs(ave_gyr.x) << "  " << devstd_gyr.y / fabs(ave_gyr.y) << "  " << devstd_gyr.z / fabs(ave_gyr.z) << std::endl;
  data_out.close();

  // Gnuplot script
  std::string gnuplot_filename = input_file.substr(0, input_file.size() - 4) + ".plt";
  std::string plot_filename = input_file.substr(0, input_file.size() - 4) + ".png";
  std::string escaped_filename = boost::replace_all_copy(input_file, "_", "\\_");

  std::ofstream plot(gnuplot_filename);
  plot << R"(reset
set terminal pngcairo dashed size 1500, 700 enhanced font 'Verdana,10'
set output ")" << plot_filename << R"("
set multiplot layout 3, 2 title 'Inertial Data Noise Test : )" << escaped_filename << R"(' font ",14"
# Styles
linew = 1.2
set style line  21 lc rgb '#0072bd' lt 7 lw linew  # blue
set style line  22 lc rgb '#d95319' lt 7 lw linew  # orange
set style line  23 lc rgb '#77ac30' lt 7 lw linew  # green
set style line  24 lc rgb '#a2142f' lt 7 lw linew  # red
set style line 101 lc rgb '#000000' lt 1 lw 1                     # black
set style line 102 lc rgb '#d6d7d9' lt 1 lw 1                     # gray
# Border xy
set border 3 front ls 101
set tics nomirror out scale 0.75
set format '%g'
set border linewidth 1.5
# Grid
set x2tics scale 0 format " "
set y2tics scale 0 format " "
set grid x2tics y2tics back ls 102
# Plot
input = ")" << input_file << R"("
set xlabel 't (s)'
set ylabel 'a_x (g)'
plot input using 10:3 every 15 notitle with lines ls 21
set xlabel 't (s)'
set ylabel '{/Symbol w}_x (dps)'
plot input using 10 : 6 every 15 notitle with lines ls 21
set xlabel 't (s)'
set ylabel 'a_y (g)'
plot input using 10 : 4 every 15 notitle with lines ls 22
set xlabel 't (s)'
set ylabel '{/Symbol w}_y (dps)'
plot input using 10 : 7 every 15 notitle with lines ls 22
set xlabel 't (s)'
set ylabel 'a_z (g)'
plot input using 10 : 5 every 15 notitle with lines ls 23
set xlabel 't (s)'
set ylabel '{/Symbol w}_z (dps)'
plot input using 10 : 8 every 15 notitle with lines ls 23
unset multiplot
)" << std::endl;
  plot.close();

  return 0;
}
