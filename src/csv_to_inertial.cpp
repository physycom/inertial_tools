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
#define _CRT_SECURE_NO_WARNINGS  // MVS C warnings shut upper
#define _SCL_SECURE_NO_WARNINGS  // MVS C warnings shut upper
#endif

#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits> 

#include "io_lib.hpp"
#include "params.h"

using namespace std;

#define MAJOR_VERSION       0
#define MINOR_VERSION       0

void usage(char * progname) {
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " -i input.csv -o output.txt" << std::endl;
}

int main(int argc, char **argv) {
  std::cout << "CSV to Inertial v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  return 0;
}