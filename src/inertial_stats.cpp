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

#include "jsoncons/json.hpp"
#include "io_lib.hpp"
#include "math_lib.h"


#define MAJOR_VERSION       0
#define MINOR_VERSION       3

void usage(char * progname) {
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cerr << "Usage: " << tokens.back() << " path/to/data/file" << std::endl;
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

  // Data parsing, conversion, storage
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
    quad_acc.yx += ay[i] * ax[i]; quad_acc.yy += ay[i] * ay[i]; quad_acc.yz += ay[i] * az[i];
    quad_acc.zx += az[i] * ax[i]; quad_acc.zy += az[i] * ay[i]; quad_acc.zz += az[i] * az[i];

    quad_gyr.xx += gx[i] * gx[i]; quad_gyr.xy += gx[i] * gy[i]; quad_gyr.xz += gx[i] * gz[i];
    quad_gyr.yx += gy[i] * gx[i]; quad_gyr.yy += gy[i] * gy[i]; quad_gyr.yz += gy[i] * gz[i];
    quad_gyr.zx += gz[i] * gx[i]; quad_gyr.zy += gz[i] * gy[i]; quad_gyr.zz += gz[i] * gz[i];

    quad_mix.xx += ax[i] * gx[i]; quad_mix.xy += ax[i] * gy[i]; quad_mix.xz += ax[i] * gz[i];
    quad_mix.yx += ay[i] * gx[i]; quad_mix.yy += ay[i] * gy[i]; quad_mix.yz += ay[i] * gz[i];
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
  cov_acc.yx = quad_acc.yx - ave_acc.y*ave_acc.x; cov_acc.yy = quad_acc.yy - ave_acc.y*ave_acc.y; cov_acc.yz = quad_acc.yz - ave_acc.y*ave_acc.z;
  cov_acc.zx = quad_acc.zx - ave_acc.z*ave_acc.x; cov_acc.zy = quad_acc.zy - ave_acc.z*ave_acc.y; cov_acc.zz = quad_acc.zz - ave_acc.z*ave_acc.z;

  cov_gyr.xx = quad_gyr.xx - ave_gyr.x*ave_gyr.x; cov_gyr.xy = quad_gyr.xy - ave_gyr.x*ave_gyr.y; cov_gyr.xz = quad_gyr.xz - ave_gyr.x*ave_gyr.z;
  cov_gyr.yx = quad_gyr.yx - ave_gyr.y*ave_gyr.x; cov_gyr.yy = quad_gyr.yy - ave_gyr.y*ave_gyr.y; cov_gyr.yz = quad_gyr.yz - ave_gyr.y*ave_gyr.z;
  cov_gyr.zx = quad_gyr.zx - ave_gyr.z*ave_gyr.x; cov_gyr.zy = quad_gyr.zy - ave_gyr.z*ave_gyr.y; cov_gyr.zz = quad_gyr.zz - ave_gyr.z*ave_gyr.z;

  cov_mix.xx = quad_mix.xx - ave_acc.x*ave_gyr.x; cov_mix.xy = quad_mix.xy - ave_acc.x*ave_gyr.y; cov_mix.xz = quad_mix.xz - ave_acc.x*ave_gyr.z;
  cov_mix.yx = quad_mix.yx - ave_acc.y*ave_gyr.x; cov_mix.yy = quad_mix.yy - ave_acc.y*ave_gyr.y; cov_mix.yz = quad_mix.yz - ave_acc.y*ave_gyr.z;
  cov_mix.zx = quad_mix.zx - ave_acc.z*ave_gyr.x; cov_mix.zy = quad_mix.zy - ave_acc.z*ave_gyr.y; cov_mix.zz = quad_mix.zz - ave_acc.z*ave_gyr.z;

  // Jsonize results
  jsoncons::json j_ave_acc;
  j_ave_acc["x"] = ave_acc.x; j_ave_acc["y"] = ave_acc.y; j_ave_acc["z"] = ave_acc.z;
  jsoncons::json j_dev_acc;
  j_dev_acc["x"] = devstd_acc.x; j_dev_acc["y"] = devstd_acc.y; j_dev_acc["z"] = devstd_acc.z;
  jsoncons::json j_quad_acc, j_quad_accx, j_quad_accy, j_quad_accz;
  j_quad_accx["x"] = quad_acc.xx; j_quad_accx["y"] = quad_acc.xy; j_quad_accx["z"] = quad_acc.xz;
  j_quad_accy["x"] = quad_acc.yx; j_quad_accy["y"] = quad_acc.yy; j_quad_accy["z"] = quad_acc.yz;
  j_quad_accz["x"] = quad_acc.zx; j_quad_accz["y"] = quad_acc.zy; j_quad_accz["z"] = quad_acc.zz;
  j_quad_acc["x"] = j_quad_accx; j_quad_acc["y"] = j_quad_accy; j_quad_acc["z"] = j_quad_accz;
  jsoncons::json j_cov_acc, j_cov_accx, j_cov_accy, j_cov_accz;
  j_cov_accx["x"] = cov_acc.xx; j_cov_accx["y"] = cov_acc.xy; j_cov_accx["z"] = cov_acc.xz;
  j_cov_accy["x"] = cov_acc.yx; j_cov_accy["y"] = cov_acc.yy; j_cov_accy["z"] = cov_acc.yz;
  j_cov_accz["x"] = cov_acc.zx; j_cov_accz["y"] = cov_acc.zy; j_cov_accz["z"] = cov_acc.zz;
  j_cov_acc["x"] = j_cov_accx; j_cov_acc["y"] = j_cov_accy; j_cov_acc["z"] = j_cov_accz;

  jsoncons::json j_ave_gyr;
  j_ave_gyr["x"] = ave_gyr.x; j_ave_gyr["y"] = ave_gyr.y; j_ave_gyr["z"] = ave_gyr.z;
  jsoncons::json j_dev_gyr;
  j_dev_gyr["x"] = devstd_gyr.x; j_dev_gyr["y"] = devstd_gyr.y; j_dev_gyr["z"] = devstd_gyr.z;
  jsoncons::json j_quad_gyr, j_quad_gyrx, j_quad_gyry, j_quad_gyrz;
  j_quad_gyrx["x"] = quad_gyr.xx; j_quad_gyrx["y"] = quad_gyr.xy; j_quad_gyrx["z"] = quad_gyr.xz;
  j_quad_gyry["x"] = quad_gyr.yx; j_quad_gyry["y"] = quad_gyr.yy; j_quad_gyry["z"] = quad_gyr.yz;
  j_quad_gyrz["x"] = quad_gyr.zx; j_quad_gyrz["y"] = quad_gyr.zy; j_quad_gyrz["z"] = quad_gyr.zz;
  j_quad_gyr["x"] = j_quad_gyrx; j_quad_gyr["y"] = j_quad_gyry; j_quad_gyr["z"] = j_quad_gyrz;
  jsoncons::json j_cov_gyr, j_cov_gyrx, j_cov_gyry, j_cov_gyrz;
  j_cov_gyrx["x"] = cov_gyr.xx; j_cov_gyrx["y"] = cov_gyr.xy; j_cov_gyrx["z"] = cov_gyr.xz;
  j_cov_gyry["x"] = cov_gyr.yx; j_cov_gyry["y"] = cov_gyr.yy; j_cov_gyry["z"] = cov_gyr.yz;
  j_cov_gyrz["x"] = cov_gyr.zx; j_cov_gyrz["y"] = cov_gyr.zy; j_cov_gyrz["z"] = cov_gyr.zz;
  j_cov_gyr["x"] = j_cov_gyrx; j_cov_gyr["y"] = j_cov_gyry; j_cov_gyr["z"] = j_cov_gyrz;

  jsoncons::json j_quad_mix, j_quad_mixx, j_quad_mixy, j_quad_mixz;
  j_quad_mixx["x"] = quad_mix.xx; j_quad_mixx["y"] = quad_mix.xy; j_quad_mixx["z"] = quad_mix.xz;
  j_quad_mixy["x"] = quad_mix.yx; j_quad_mixy["y"] = quad_mix.yy; j_quad_mixy["z"] = quad_mix.yz;
  j_quad_mixz["x"] = quad_mix.zx; j_quad_mixz["y"] = quad_mix.zy; j_quad_mixz["z"] = quad_mix.zz;
  j_quad_mix["x"] = j_quad_mixx; j_quad_mix["y"] = j_quad_mixy; j_quad_mix["z"] = j_quad_mixz;
  jsoncons::json j_cov_mix, j_cov_mixx, j_cov_mixy, j_cov_mixz;
  j_cov_mixx["x"] = cov_mix.xx; j_cov_mixx["y"] = cov_mix.xy; j_cov_mixx["z"] = cov_mix.xz;
  j_cov_mixy["x"] = cov_mix.yx; j_cov_mixy["y"] = cov_mix.yy; j_cov_mixy["z"] = cov_mix.yz;
  j_cov_mixz["x"] = cov_mix.zx; j_cov_mixz["y"] = cov_mix.zy; j_cov_mixz["z"] = cov_mix.zz;
  j_cov_mix["x"] = j_cov_mixx; j_cov_mix["y"] = j_cov_mixy; j_cov_mix["z"] = j_cov_mixz;

  jsoncons::json jresults;
  jresults["data_samples"] = ax.size();
  jresults["acc_ave_g"] = j_ave_acc;
  jresults["acc_devstd_g"] = j_dev_acc;
  jresults["acc_quad_g^2"] = j_quad_acc;
  jresults["acc_cov_g^2"] = j_cov_acc;
  jresults["gyr_ave_g"] = j_ave_gyr;
  jresults["gyr_devstd_g"] = j_dev_gyr;
  jresults["gyr_quad_g^2"] = j_quad_gyr;
  jresults["gyr_cov_g^2"] = j_cov_gyr;
  jresults["mix_quad_g^2"] = j_quad_mix;
  jresults["mix_cov_g^2"] = j_cov_mix;

  // Output
  std::string results_filename = input_file.substr(0, input_file.size() - 4) + "_stats.json";
  std::string gnuplot_filename = input_file.substr(0, input_file.size() - 4) + "_stats.plt";
  std::string plot_filename = input_file.substr(0, input_file.size() - 4) + "_stats.png";
  std::string escaped_filename = boost::replace_all_copy(input_file, "_", "\\_");

  std::ofstream results(results_filename);
  results << jsoncons::pretty_print(jresults) << std::endl;
  results.close();

  std::ofstream plot(gnuplot_filename);
  plot << R"(reset
set terminal pngcairo dashed size 1500, 700 enhanced font 'Verdana,10'
set output ")" << plot_filename << R"("
set multiplot layout 3, 2 title 'Inertial Data Noise Test : )" << escaped_filename << R"(' font ",14"
set key opaque
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
# acc_x
set xlabel 't (s)'
set ylabel 'a_x (g)'
plot input using 10:3 every 15 title '{/Symbol s}_{a_x} = )" << fixed << setprecision(3) << devstd_acc.x << R"( (g)' with lines ls 21
# gyr_x
set xlabel 't (s)'
set ylabel '{/Symbol w}_x (dps)'
plot input using 10:6 every 15 title '{/Symbol s}_{{/Symbol w}_x} = )" << fixed << setprecision(2) << devstd_gyr.x << R"( (dps)' with lines ls 21
# acc_y
set xlabel 't (s)'
set ylabel 'a_y (g)'
plot input using 10:4 every 15 title '{/Symbol s}_{a_y} = )" << fixed << setprecision(3) << devstd_acc.y << R"( (g)' with lines ls 22
# gyr_y
set xlabel 't (s)'
set ylabel '{/Symbol w}_y (dps)'
plot input using 10:7 every 15 title '{/Symbol s}_{{/Symbol w}_y} = )" << fixed << setprecision(2) << devstd_gyr.y << R"( (dps)' with lines ls 22
# acc_z
set xlabel 't (s)'
set ylabel 'a_z (g)'
plot input using 10:5 every 15 title '{/Symbol s}_{a_z} = )" << fixed << setprecision(3) << devstd_acc.z << R"( (g)' with lines ls 23
# gyr_z
set xlabel 't (s)'
set ylabel '{/Symbol w}_z (dps)'
plot input using 10:8 every 15 title '{/Symbol s}_{{/Symbol w}_z} = )" << fixed << setprecision(2) << devstd_gyr.z << R"( (dps)' with lines ls 23
unset multiplot
)" << std::endl;
  plot.close();

  return 0;
}
