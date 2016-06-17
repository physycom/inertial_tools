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

#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits> 

#include "io_lib.hpp"
#include "params.h"
#include "math_func.h"

using namespace std;

#define MAJOR_VERSION       1
#define MINOR_VERSION       0

void usage(char * progname) {
  std::cout << "Usage: " << progname << " path/to/data/file" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file, compliant with PHYSYCOM inertial standard" << std::endl;
  exit(-3);
}

int main(int argc, char **argv) {
  std::cout << "Speed Compare v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  std::string input_file;
  if (argc == 2) {
    input_file = argv[1];
    if (input_file.substr(0, 2) == ".\\") input_file = input_file.substr(2, input_file.size() - 2);
  }
  else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
  }

  vector< vector<string> > file_tokens;
  file_tokens = Read_from_file(input_file);
  vector< vector<double> > data = tokens_to_double(file_tokens);

  vector<double> gps_speed, inertial_speed, times, ax, v;
  for (size_t i = 0; i < data.size(); i++) {
    times.push_back(data[i][TIMESTAMP_INDEX]);
    gps_speed.push_back(data[i][SPEED_INDEX]);
    ax.push_back(data[i][AX_INDEX] * 9.81);
  }

  inertial_speed = Integrate(times, ax, 0);

  std::ofstream output("compare.txt");
  for (size_t i = 0; i < inertial_speed.size(); i++) {
    output << std::fixed << std::setprecision(4)
      << times[i] - times[0] << "\t" << inertial_speed[i]/3.6 << "\t" << gps_speed[i] << std::endl;
  }
  output.close();

  return 0;
}