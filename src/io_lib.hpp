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

#include <boost/algorithm/string.hpp>
#include <boost/utility.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

bool Belongs_to_string(char c, std::string s)
{
  for (size_t i = 0; i < s.size(); i++) {
    if (c == s.at(i))
      return true;
  }
  return false;
}

std::vector<std::vector<std::string>> Read_from_file(std::string file_name)
{
  std::string line;
  std::vector<std::string> tokens;
  std::vector<std::vector<std::string>> file_tokens;
  std::ifstream data_file(file_name);
  if (!data_file) {
    std::cerr << "ERROR: file " << file_name << " could not be opened" << std::endl;
    exit(-1);
  }
  while (!data_file.eof()) {
    line.clear();
    tokens.clear();
    std::getline(data_file, line);
    if (line.size() > 0) {
      boost::algorithm::split(tokens, line, boost::algorithm::is_any_of("\t"));
      if (tokens[0].at(0) == '#')
        continue;
      else
        file_tokens.push_back(tokens);
    }
  }
  return file_tokens;
}

std::vector<std::vector<double>> tokens_to_double(std::vector<std::vector<std::string>> parsed_file)
{
  std::vector<double> doubled_line;
  std::vector<std::vector<double>> doubled_file;

  for (auto& i : parsed_file) {
    bool skip_line = false;
    doubled_line.clear();
    doubled_line.resize(i.size());
    for (size_t j = 0; j < i.size(); j++) {
      try {
        doubled_line[j] = stod(i[j]);
      } catch (...) { // to avoid line containing empty values
        skip_line = true;
        break;
      }
      if (std::isnan(doubled_line[j])) { // to avoid line containing NaN values
        skip_line = true;
        break;
      }
    }
    if (!skip_line)
      doubled_file.push_back(doubled_line);
  }
  return doubled_file;
}

void dump_to_csv(std::vector<std::vector<double>> data, std::string filename)
{
  std::ofstream output(filename);
  if (!output) {
    std::cerr << "FAILED: file " << filename << " could not be opened" << std::endl;
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
