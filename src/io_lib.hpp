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



#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

#include <boost/utility.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost::algorithm;

bool Belongs_to_string(char c, string s) {
  for (size_t i = 0; i < s.size(); i++) {
    if (c == s.at(i)) return true;
  }
  return false;
}

vector< vector<string> > Read_from_file(string file_name) {
  string line; vector<string> tokens; vector< vector<string> > file_tokens;
  fstream data_file; data_file.open(file_name.c_str());
  if (!data_file) {
    cout << "FAILED: file " << file_name << " could not be opened" << endl << "Press enter to close.";
    cin.get();
    exit(777);
  }
  while (!data_file.eof()) {
    line.clear();
    tokens.clear();
    getline(data_file, line);
    if (line.size() > 0) {
//      inertial files only have \t separation
//      split(tokens, line, is_any_of("\t ,;:"));
      split(tokens, line, is_any_of("\t"));
      if (tokens[0].at(0) == '#') continue;
      else file_tokens.push_back(tokens);
    }
  }
  return file_tokens;
}

std::vector< std::vector<double> > tokens_to_double(std::vector< std::vector<std::string> > parsed_file) {
  std::vector<double> doubled_line;
  std::vector< std::vector<double> > doubled_file;

  for (auto &i : parsed_file) {
    doubled_line.clear();
    doubled_line.resize(i.size());
    for (size_t j = 0; j < i.size(); j++) doubled_line[j] = atof(i[j].c_str());
    doubled_file.push_back(doubled_line);
  }
  return doubled_file;
}


void dump_to_csv(std::vector< std::vector<double> > data, std::string filename) {
  std::ofstream output(filename);
  if (!output) {
    std::cout << "FAILED: file " << filename << " could not be opened" << std::endl << "Press enter to close.";
    std::cin.get();
    exit(-2);
  }

  for (size_t i = 0; i < data.size(); i++) {
    for (size_t j = 0; j < data[i].size(); j++) {
      output << std::fixed << std::setprecision(3) << std::setw(7) << data[i][j] << '\t';
    }
    output << std::endl;
  }
  output.close();

  return;
}
