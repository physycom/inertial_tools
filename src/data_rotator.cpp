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
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "io_lib.hpp"
#include "libbbmutils/bbmutils.h"
#include "math_func.h"
#include "params.h"

#define MAJOR_VERSION 2
#define MINOR_VERSION 0

#define HORIZONTAL_MODE 'h'
#define VERTICAL_MODE 'v'
#define BOTH_MODE 'b'

void usage(char* progname)
{
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " alpha nx ny nz path/to/data/file" << std::endl;
  std::cout << "       alpha : angle of rotation in DEGREES" << std::endl;
  std::cout << "       (nx, ny, nz) : components of the rotation axis" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file, compliant with PHYSYCOM inertial standard" << std::endl;
  exit(-3);
}

int main(int argc, char** argv)
{
  std::cout << "Data Rotator v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl
            << std::endl;

  std::string input_file;
  double alpha, nx, ny, nz;
  if (argc == 6) {
    input_file = argv[5];
    try {
      alpha = std::stod(std::string(argv[1])) / RAD_TO_DEG;
      nx = std::stod(std::string(argv[2]));
      ny = std::stod(std::string(argv[3]));
      nz = std::stod(std::string(argv[4]));
    } catch (std::exception& e) {
      std::cout << "EXCEPTION: " << e.what() << std::endl;
      usage(argv[0]);
    }
  } else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
  }

  std::cout << "Rotating file : " << input_file << "\tAngle : " << alpha * RAD_TO_DEG
            << "\tAxis : ( " << nx << " , " << ny << " , " << nz << " ) " << std::endl;

  std::vector<std::vector<std::string>> file_tokens = Read_from_file(input_file);
  std::vector<std::vector<double>> data = tokens_to_double(file_tokens);

  VEC3D axis;
  set_vec3d(&axis, nx, ny, nz);
  MAT3D rotation;
  make_rotation(&rotation, &axis, alpha);
  std::vector<std::vector<double>> data_r = rotate_inertial(data, rotation);
  std::string outfile = input_file.substr(0, input_file.size() - 4) + "_rot.txt";
  dump_to_csv(data_r, outfile);

  return 0;
}
