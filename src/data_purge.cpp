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

#include "io_lib.hpp"

using namespace std;

#define MAJOR_VERSION       3
#define MINOR_VERSION       0


//////// Filter object
class Filter {
public:
  size_t first, last;
  size_t left, right;
  std::vector<double> params;

  Filter(const size_t& data_size, const std::string& filter_path) {
    std::vector<double> coeff;
    std::ifstream filter_file(filter_path);
    if (!filter_file) {
      cerr << "ERROR: file " << filter_path << " could not be opened" << endl;
      exit(1);
    }
    while (!filter_file.eof()) {
      double temp;
      filter_file >> temp;
      coeff.push_back(temp);
    }
    filter_file.close();

    double normalization = 0.0;
    for (auto p : coeff) normalization += p;
    for (auto p : coeff) params.push_back(p / normalization);

    first = (params.size() % 2) ? (params.size() / 2) : (params.size() / 2 - 1);
    last = data_size - params.size() / 2 - 1;
    left = (params.size() % 2) ? (params.size() / 2) : (params.size() / 2 - 1);
    right = params.size() / 2;
  };
};

std::ostream& operator<<(std::ostream &ost, const Filter &f) {
  ost << "first = " << f.first << "  last = " << f.last << "  size = " << f.params.size() << "  left = " << f.left << "  right = " << f.right;
  return ost;
}


//////// MAIN
void usage(char * progname) {
  std::vector<std::string> tokens;
  boost::split(tokens, std::string(progname), boost::is_any_of("/\\"));
  std::cout << "Usage: " << tokens.back() << " -filter path/to/filter [-subtract] path/to/data/file" << std::endl;
  std::cout << "       -filter create a filter from file" << std::endl;
  std::cout << "       -subtract remove the average value of each column from the data [optional]" << std::endl;
  std::cout << "       path/to/data/file must be a tab separated value file with timestamps in the first column" << std::endl;
}

int main(int argc, char **argv) {
  std::cout << "Data Purge v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  std::string input_file, filter_file;
  int ma_sample = 10;
  bool subtract_ave = false;
  if (argc == 4 || argc == 5) {
    for (int i = 1; i < argc; i++) {
      if (std::string(argv[i]) == "-filter") filter_file = argv[++i];
      else if (std::string(argv[i]) == "-subtract") subtract_ave = true;
      else input_file = argv[i];
    }
  }
  else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
    exit(-1);
  }

  if (input_file.substr(0, 2) == ".\\" || input_file.substr(0, 2) == "./") input_file = input_file.substr(2, input_file.size() - 2);
  std::cout << "Purging        : " << input_file << std::endl;
  std::cout << "Filter         : " << filter_file << std::endl;
  std::cout << "Remove average : " << (subtract_ave ? "ON" : "OFF") << std::endl;

  // Parse text file and convert it to vector of double
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

  // Create filter object from file
  Filter f(data.size(), filter_file);

  // Remove global offset
  if (subtract_ave) {
    vector<double> averages(data[0].size(), 0.0);
    // Calculate averages for each column, except the first
    for (size_t i = 1; i < data[0].size(); i++) {
      for (size_t j = 0; j < data.size(); j++) {
        averages[i] += data[j][i];
      }
    }
    // Normalize it
    for (size_t i = 0; i < averages.size(); i++) {
      averages[i] /= (double)data.size();
    }
    // Subtract from the whole set of data
    for (size_t i = 0; i < data.size(); i++) {
      for (size_t j = 0; j < data[i].size(); j++) {
        data[i][j] -= averages[j];
      }
    }
  }

  // Filtering
  for (auto p : f.params) std::cout << p << " "; std::cout << std::endl;

  vector< vector<double> > data_purged;
  for (size_t i = f.first; i <= f.last; i++) {
    vector<double> row_purged;
    row_purged.push_back(data[i][0]);                   // pushback timestamp without modification
    for (size_t j = 1; j < data[0].size(); j++) {
      double value_purged = 0.0;
      for (size_t k = i - f.left; k < i + f.right; k++) {
        value_purged += data[k][j] * f.params[k - i + f.left];       // params[...] the index is calculated to keep it in range
      }
      row_purged.push_back(value_purged);
    }
    data_purged.push_back(row_purged);
  }
  std::cout << "Data original : r " << data.size() << " c " << data[0].size() << std::endl;
  std::cout << "Data purged   : r " << data_purged.size() << " c " << data_purged[0].size() << std::endl;


  // Dump to file
  std::string out_file = input_file.substr(0, input_file.size() - 4) + "_purged" + (subtract_ave ? "_sub" : "") + ".txt";
  std::ofstream output(out_file);
  output << "# " << f << std::endl;
  for (size_t i = 0; i < data_purged.size(); i++) {
    for (size_t j = 0; j < data_purged[i].size(); j++) {
      output << std::fixed << std::setprecision(6) << data_purged[i][j] << "\t";
    }
    output << std::endl;
  }
  output.close();
  std::cout << "Purge result in : " << out_file << std::endl;

  return 0;
}