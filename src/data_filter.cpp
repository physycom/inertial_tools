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

#define MAJOR_VERSION       2
#define MINOR_VERSION       0

void usage(char * progname) {
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " path/to/data/file column_index min_value max_value" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file" << std::endl;
  std::cout << "       column_index is an int indicating the column to filter, starting from 0" << std::endl;
  std::cout << "       min_value (max_value) represents the min (max) value which passes the filtering, use \"unset\" for both to unset" << std::endl;
}

int main(int argc, char **argv) {
  std::cout << "Data Filter v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  std::string input_file;
  int col_index;
  double min, max;
  if (argc > 4) {
    input_file = argv[1];
    if (input_file.substr(0, 2) == ".\\") input_file = input_file.substr(2, input_file.size() - 2);

    try {
      col_index = std::stoi(std::string(argv[2]));
    }
    catch (std::exception &e) {
      std::cout << "Exception thrown in parsing column index : " << e.what() << std::endl;
      exit(-1);
    }

    try {
      min = std::stod(std::string(argv[3]));
    }
    catch (std::exception &e) {
      if (std::string(argv[3]) == "unset") min = -std::numeric_limits<double>::max();
      else {
        std::cout << "Exception thrown in parsing min value : " << e.what() << std::endl;
        exit(-1);
      }
    }

    try {
      max = std::stod(std::string(argv[4]));
    }
    catch (std::exception &e) {
      if (std::string(argv[4]) == "unset") max = std::numeric_limits<double>::max();
      else {
        std::cout << "Exception thrown in parsing max value : " << e.what() << std::endl;
        exit(-1);
      }
    }
    
  }
  else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
    exit(-3);
  }

  std::cout << "Filtering file : " << input_file << "\tColumn : " << col_index 
    << "\tRange : [ " << min << " , " << max << " ]" << std::endl;

  vector< vector<string> > file_tokens;
  file_tokens = Read_from_file(input_file);
  vector< vector<double> > data;
  for (size_t i = 0; i < file_tokens.size(); i++) {
    vector<double> row;
    for (size_t j = 0; j < file_tokens[i].size(); j++) {
      row.push_back(std::stod(file_tokens[i][j]));
    }
    data.push_back(row);
  }

  vector< vector<double> > data_filtered;
  for (size_t i = 0; i < data.size(); i++) {
    if (data[i][col_index] > min && data[i][col_index] < max) data_filtered.push_back(data[i]);
  }

  std::cout << "Filtering complete - original lines : " << data.size() << "\t filtered lines : " << data_filtered.size() 
    << " (" << int(100*data_filtered.size()/double(data.size())) << " %)" << std::endl;

  std::string output = input_file.substr(0, input_file.size() - 4) + "_filtered.txt";
  dump_to_csv(data_filtered, output);
  std::cout << "Filter result dumped to : \"" << output << "\"" << std::endl;

  return 0;
}