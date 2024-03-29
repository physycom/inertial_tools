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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "io_lib.hpp"
#include <bbmutils.h>
#include "math_func.h"
#include "params.h"

#include <boost/algorithm/string.hpp>
#include <boost/utility.hpp>

#define MAJOR_VERSION 1
#define MINOR_VERSION 2

#define HORIZONTAL_MODE 'h'
#define VERTICAL_MODE 'v'
#define BOTH_MODE 'b'

#define V_MIN 10.0 /* km/h */
#define OMEGA_MIN 5.0 /* dps  */
#define OMEGA_MAX 10.0 /* dps  */
#define A_MAX 1.5 /* g    */

//////// CALIBRATION ALGORITHM
void find_angle_V(double* thetaV, double* dthetaV, VEC3D* axis, std::vector<std::vector<double>> data)
{
  int V_samples = 0;
  VEC3D accV_mean;
  memset(&accV_mean, 0, sizeof(VEC3D));
  MAT3D accV_cov;
  memset(&accV_cov, 0, sizeof(MAT3D));

  // filter data
  for (size_t i = 0; i < data.size(); i++) {
    double ax = data[i][AX_INDEX],
           ay = data[i][AY_INDEX],
           az = data[i][AZ_INDEX],
           gx = data[i][GX_INDEX],
           gy = data[i][GY_INDEX],
           gz = data[i][GZ_INDEX];

    double a_squared = ax * ax + ay * ay + az * az;
    double g_squared = gx * gx + gy * gy + gz * gz;

    if (a_squared < A_MAX * A_MAX && g_squared < OMEGA_MAX * OMEGA_MAX) /* stationary regime  */
    {
      /* data stored for vertical angle */
      accV_mean.x += (float)ax;
      accV_mean.y += (float)ay;
      accV_mean.z += (float)az;
      accV_cov.xx += float(ax * ax);
      accV_cov.xy += float(ax * ay);
      accV_cov.xz += float(ax * az);
      accV_cov.yx += float(ay * ax);
      accV_cov.yy += float(ay * ay);
      accV_cov.yz += float(ay * az);
      accV_cov.zx += float(az * ax);
      accV_cov.zy += float(az * ay);
      accV_cov.zz += float(az * az);
      V_samples++;
    }
  }
  multiply_vec3d(&accV_mean, 1.0 / (float)V_samples);
  multiply_mat3d(&accV_cov, 1.0 / (float)V_samples);

  std::cout << "Evaluating VERTICAL angle from " << V_samples << "/" << data.size() << " ( " << int(100 * V_samples / double(data.size())) << " %) samples" << std::endl;

  /* axis */
  VEC3D ortho, acc_n, z_axis;
  set_vec3d(&z_axis, 0., 0., 1.);
  memcpy(&acc_n, &accV_mean, sizeof accV_mean);
  normalize_vec3d(&acc_n);
  prod_cross_3d(&ortho, &acc_n, &z_axis);
  memcpy(&acc_n, &accV_mean, sizeof accV_mean);
  memcpy(axis, &ortho, sizeof ortho);
  normalize_vec3d(axis);

  /* angle */
  double costhetaV, sinthetaV, errthetaV;
  costhetaV = prod_dot_3d(&acc_n, &z_axis);
  sinthetaV = ortho.mod;
  errthetaV = 1 - costhetaV * costhetaV - sinthetaV * sinthetaV;
  *thetaV = atan2(sinthetaV, costhetaV);

  /* error */
  double x, y, z, dx, dy, dz;
  x = accV_mean.x;
  y = accV_mean.y;
  z = accV_mean.z;
  dx = sqrt(accV_cov.xx - x * x);
  dy = sqrt(accV_cov.yy - y * y);
  dz = sqrt(accV_cov.zz - z * z);

  double c, dc, s, ds;

  c = costhetaV;
  s = sinthetaV;

  dc = (z * x * dx) * (z * x * dx) + (z * y * dy) * (z * y * dy) + (x * x + y * y) * (x * x + y * y) * dz * dz;
  dc = sqrt(dc);
  dc /= (sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z) * sqrt(x * x + y * y + z * z));

  ds = (z * z * x * dx) * (z * z * x * dx) / ((x * x + y * y) * (x * x + y * y)) + (z * z * y * dy) * (z * z * y * dy) / ((x * x + y * y) * (x * x + y * y)) + z * z * dz * dz;
  ds = sqrt(ds);
  ds *= sqrt((x * x + y * y) / ((x * x + y * y + z * z) * (x * x + y * y + z * z) * (x * x + y * y + z * z)));

  *dthetaV = s * s * dc * dc + c * c * ds * ds;
  *dthetaV = sqrt(*dthetaV);
  *dthetaV /= (s * s + c * c);

  std::cout << "Vertical angle (deg) " << *thetaV * RAD_TO_DEG << " +- " << *dthetaV * RAD_TO_DEG << std::endl;
  std::cout << "Axis ( " << axis->x << " , " << axis->y << " , " << axis->z << " ) " << std::endl;

  return;
}

void find_angle_H(double* thetaH, double* dthetaH, std::vector<std::vector<double>> data)
{
  int H_samples = 0;
  VEC3D accH_mean;
  memset(&accH_mean, 0, sizeof(VEC3D));
  MAT3D accH_cov;
  memset(&accH_cov, 0, sizeof(MAT3D));

  // filter data
  for (size_t i = 0; i < data.size(); i++) {
    double v = data[i][SPEED_INDEX],
           ax = data[i][AX_INDEX],
           ay = data[i][AY_INDEX],
           az = data[i][AZ_INDEX],
           gx = data[i][GX_INDEX],
           gy = data[i][GY_INDEX],
           gz = data[i][GZ_INDEX];

    if (v > V_MIN && gz > OMEGA_MIN) /* fast driving regime */
    {
      double sign = (gz > 0) - (gz < 0);
      ax = sign * ax;
      ay = sign * ay;
      az = az;

      /* data stored for horizontal angle  */
      accH_mean.x += ax;
      accH_mean.y += ay;
      accH_mean.z += az;
      accH_cov.xx += ax * ax;
      accH_cov.xy += ax * ay;
      accH_cov.xz += ax * az;
      accH_cov.yx += ay * ax;
      accH_cov.yy += ay * ay;
      accH_cov.yz += ay * az;
      accH_cov.zx += az * ax;
      accH_cov.zy += az * ay;
      accH_cov.zz += az * az;

      H_samples++;
    }
  }
  multiply_vec3d(&accH_mean, 1.0 / (float)H_samples);
  multiply_mat3d(&accH_cov, 1.0 / (float)H_samples);

  std::cout << "Evaluating HORIZONTAL angle from " << H_samples << "/" << data.size() << " ( " << int(100 * H_samples / double(data.size())) << " %) samples" << std::endl;

  /* angle */
  double costhetaH, sinthetaH, errthetaH, x, y, dx, dy;
  x = accH_mean.x;
  y = accH_mean.y;
  costhetaH = y / sqrt(x * x + y * y); /* + < v , y > */
  sinthetaH = x / sqrt(x * x + y * y); /* - < v , x > */
  errthetaH = 1 - costhetaH * costhetaH - sinthetaH * sinthetaH;
  *thetaH = atan2(costhetaH, sinthetaH);
  if (*thetaH < 0)
    *thetaH += 2 * M_PI;

  /* error */
  double nx, dnx, ny, dny, axy, d_axy, r, dr;
  dx = sqrt(accH_cov.xx - x * x);
  dy = sqrt(accH_cov.yy - y * y);
  x = fabs(x);
  y = fabs(y);
  axy = sqrt(x * x + y * y);
  d_axy = (x * dx + y * dy) / axy;
  nx = x / axy;
  dnx = dx / axy + x / axy / axy * d_axy;
  ny = y / axy;
  dny = dy / axy + y / axy / axy * d_axy;
  r = nx / ny;
  dr = dnx / ny + nx / ny / ny * dny;
  *dthetaH = dr / (1 + r * r);

  /* output */
  std::cout << "Horizontal angle (deg) " << *thetaH * RAD_TO_DEG << " +- " << *dthetaH * RAD_TO_DEG << std::endl;

  return;
}

void usage(char* progname)
{
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " -[h/v/b] path/to/data/file " << std::endl;
  std::cout << "       [h]orizontal calibration" << std::endl;
  std::cout << "       [v]ertical calibration" << std::endl;
  std::cout << "       [b]oth calibration" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file, compliant with PHYSYCOM inertial standard" << std::endl;
  exit(-3);
}

int main(int argc, char** argv)
{
  std::cout << "Calibrator v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl
            << std::endl;

  std::string input_file;
  char mode;
  if (argc == 3) {
    mode = argv[1][1];
    input_file = argv[2];
  } else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
  }

  std::cout << "Calibrating from file : " << input_file << "\tMode : " << mode << std::endl;

  std::vector<std::vector<std::string>> file_tokens = Read_from_file(input_file);
  std::vector<std::vector<double>> data = tokens_to_double(file_tokens);

  if (mode == HORIZONTAL_MODE) {
    double thetaH = 0, dthetaH = 0;
    VEC3D z_axis;
    set_vec3d(&z_axis, 0.0, 0.0, 1.0);
    find_angle_H(&thetaH, &dthetaH, data);
    MAT3D rotation;
    make_rotation(&rotation, &z_axis, thetaH);
    std::vector<std::vector<double>> data_r = rotate_inertial(data, rotation);
    std::string outfile = input_file.substr(0, input_file.size() - 4) + "_rotH.txt";
    dump_to_csv(data_r, outfile);
  } else if (mode == VERTICAL_MODE) {
    double thetaV = 0, dthetaV = 0;
    VEC3D axis;
    find_angle_V(&thetaV, &dthetaV, &axis, data);
    MAT3D rotation;
    make_rotation(&rotation, &axis, thetaV);
    std::vector<std::vector<double>> data_r = rotate_inertial(data, rotation);
    std::string outfile = input_file.substr(0, input_file.size() - 4) + "_rotV.txt";
    dump_to_csv(data_r, outfile);
  } else if (mode == BOTH_MODE) {
    // coming soon
  } else {
    std::cout << "Mode : " << mode << " unknown" << std::endl;
    usage(argv[0]);
  }

  return 0;
}
