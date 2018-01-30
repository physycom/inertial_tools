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
#include <string>
#include <vector>

#include "io_lib.hpp"
#include "params.h"

#define MAJOR_VERSION 0
#define MINOR_VERSION 0

void usage(char* progname)
{
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " -i input.csv -o inertial.txt" << std::endl;
}

int main(int argc, char** argv)
{
  std::cout << "CSV to Inertial v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl
            << std::endl;

  std::string input_name, output_name("");
  if (argc > 2) { /* Parse arguments, if there are arguments supplied */
    for (int i = 1; i < argc; i++) {
      if ((argv[i][0] == '-') || (argv[i][0] == '/')) { // switches or options...
        switch (tolower(argv[i][1])) {
        case 'i':
          input_name = argv[++i];
          break;
        case 'o':
          output_name = argv[++i];
          break;
        default: // no match...
          std::cerr << "ERROR: Flag \"" << argv[i] << "\" not recognized. Quitting..." << std::endl;
          usage(argv[0]);
          exit(1);
        }
      } else {
        std::cerr << "ERROR: Flag \"" << argv[i] << "\" not recognized. Quitting..." << std::endl;
        usage(argv[0]);
        exit(2);
      }
    }
  } else {
    std::cerr << "ERROR: No flags specified. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
    exit(3);
  }

  // Safety and improvements in file names
  if (input_name == "") {
    std::cerr << "ERROR: No input file specified. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
    exit(4);
  } else if (input_name.substr(0, 2) == "./")
    input_name = input_name.substr(2, input_name.size() - 2);
  if (output_name == "")
    output_name = "inertial.txt";

  // Courtesy cout
  std::cout << "Input CSV       : " << input_name << std::endl;
  std::cout << "Output inertial : " << output_name << std::endl;

  // Parsing input
  double _ax, _ay, _az;
  std::vector<std::vector<double>> acc;
  std::ifstream input(input_name);
  if (!input) {
    std::cerr << "ERROR: Unable to open " << input_name << ". Quitting..." << std::endl;
    exit(5);
  }
  while (input >> _ax >> _ay >> _az) { // NO SAFETY because mattia says so
    acc.push_back({ _ax, _ay, _az });
  }
  input.close();

  std::ofstream output(output_name);
  if (!output) {
    std::cerr << "ERROR: Unable to create " << output_name << ". Quitting..." << std::endl;
    exit(6);
  }
  output << "#  timestamp # speed #     ax #     ay #     az #     gx #      gy #     gz #   |a| #  timestamp_rel #" << std::endl;
  output << "#            #       #        #        #        #        #         #        #       #                #" << std::endl;
  double t0 = -0.01;
  for (size_t i = 0; i < acc.size(); ++i) {
    output << std::fixed << std::setprecision(4) << (t0 += 0.01) << "\t" << 0 << "\t"
           << std::setprecision(5) << acc[i][0] << "\t" << acc[i][1] << "\t" << acc[i][2] << "\t"
           << std::setprecision(5) << 0 << "\t" << 0 << "\t" << 0 << "\t"
           << std::setprecision(5) << 0 << "\t"
           << std::setprecision(4) << 0 << std::endl;
  }
  output.close();

  return 0;
}
