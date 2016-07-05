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
#include <algorithm>
#include <limits>

#include "io_lib.hpp"
#include "params.h"

using namespace std;

#define GRID_STEP           500

//////// Math functions
double entropy(vector<double> p) {
  double e = 0.;
  for (auto pi : p) if (pi > numeric_limits<double>::epsilon()) e += pi*log(pi);
  return -e;
}

double mutual_entropy(std::vector<std::vector<double>> pxy, std::vector<double> px, std::vector<double> py) {
  double me = 0.;
  for (size_t i = 0; i < pxy.size(); i++) {
    for (size_t j = 0; j < pxy[i].size(); j++) {
      me += ( (pxy[i][j] > numeric_limits<double>::epsilon()) ? pxy[i][j] * log(pxy[i][j] / (px[i] * py[j])) : 0.);
    }
  }

  return me;
}

double mutual_entropy(std::vector<std::vector<double>> pxy) {
  std::vector<double> px(pxy.size(), 0.), py(pxy[0].size(), 0.);
  for (size_t i = 0; i < pxy.size(); i++) {
    for (size_t j = 0; j < pxy[i].size(); j++) {
      px[i] += pxy[i][j];
      py[j] += pxy[i][j];
    }
  }

  return mutual_entropy(pxy, px, py);
}


//////// MAIN
#define MAJOR_VERSION       2
#define MINOR_VERSION       0

#define DYNAMIC_BIN         1
#define FIXED_BIN           2

void usage(char * progname) {
  std::cout << "Usage: " << progname << " -bf <bin_fraction> -s <index_shift> path/to/data/file" << std::endl;
  std::cout << "       DINAMIC BINNING: <bin_fraction> is a decimal representing the fraction of data in each bin (e.g. 0.05 corresponds to 5%)" << std::endl;
  std::cout << "       <index_shift> is an integer representing the shift in indices between GZ and AX,AY" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file, compliant with PHYSYCOM inertial standard" << std::endl << std::endl;

  std::cout << "Usage: " << progname << " -bn <bin_number> -s <index_shift> path/to/data/file" << std::endl;
  std::cout << "       FIXED BINNING: <bin_number> is an integer representing the number of bins for each set of data" << std::endl;
  std::cout << "       <index_shift> is an integer representing the shift in indices between GZ and AX,AY" << std::endl;
  std::cout << "       path/to/data/file must be a tab-separated value file, compliant with PHYSYCOM inertial standard" << std::endl << std::endl;

  exit(-3);
}

int main(int argc, char **argv) {
  std::cout << "Mutual Entropy Calculator v" << MAJOR_VERSION << "." << MINOR_VERSION << std::endl << std::endl;

  std::string input_file;
  double bin_fraction;
  int bin_number, index_shift = 0;
  char mode;
  if (argc > 3) {
    for (int i = 0; i < argc; i++) {
      if (std::string(argv[i]) == "-bf") {
        mode = DYNAMIC_BIN;
        try { bin_fraction = std::stod(std::string(argv[++i])); }
        catch (std::exception &e) { std::cout << "EXCEPTION: " << e.what() << std::endl; usage(argv[0]); }
      }
      if (std::string(argv[i]) == "-bn") {
        mode = FIXED_BIN;
        try { bin_number = std::stoi(std::string(argv[++i])); }
        catch (std::exception &e) { std::cout << "EXCEPTION: " << e.what() << std::endl; usage(argv[0]); }
      }
      if (std::string(argv[i]) == "-s") {
        try { index_shift = std::stoi(std::string(argv[++i])); }
        catch (std::exception &e) { std::cout << "EXCEPTION: " << e.what() << std::endl; usage(argv[0]); }
      }
    }
    input_file = argv[argc - 1];
    if (input_file.substr(0, 2) == ".\\") input_file = input_file.substr(2, input_file.size() - 2);
  }
  else {
    std::cout << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << std::endl;
    usage(argv[0]);
  }

  // data parsing, convertion, storage
  vector< vector<string> > file_tokens = Read_from_file(input_file);
  vector< vector<double> > data = tokens_to_double(file_tokens);
  vector<double> ax_, ay_, gz_, ax, ay, gz;
  for (auto line : data) {
    ax_.push_back(line[AX_INDEX]);
    ay_.push_back(line[AY_INDEX]);
    gz_.push_back(line[GZ_INDEX]);
  }

  // shift, if any
  if (index_shift > (int)ax_.size() / 2) {
    std::cout << "ERROR: Index shift too big, upper bound : " << ax_.size() / 2 << std::endl;
    exit(21);
  }
  else {
    std::cout << "SHIFT MODE (index shift " << index_shift << " ) " << std::endl;
  }
  for (size_t i = 0; i < ax_.size() - abs(index_shift); i++) {
    if (index_shift >= 0) {
      ax.push_back(ax_[i]);
      ay.push_back(ay_[i]);
      gz.push_back(gz_[i + index_shift]);
    }
    else {
      ax.push_back(ax_[i - index_shift]);
      ay.push_back(ay_[i - index_shift]);
      gz.push_back(gz_[i]);
    }
  }

  // ranges
  auto   ax_mm = std::minmax_element(ax.begin(), ax.end());
  double ax_min = *ax_mm.first;
  double ax_max = *ax_mm.second;

  auto   ay_mm = std::minmax_element(ay.begin(), ay.end());
  double ay_min = *ay_mm.first;
  double ay_max = *ay_mm.second;

  auto   gz_mm = std::minmax_element(gz.begin(), gz.end());
  double gz_min = *gz_mm.first;
  double gz_max = *gz_mm.second;

  // binning
  vector<vector<int>> counter_ax_gz, counter_ay_gz;
  switch (mode) {
  case DYNAMIC_BIN:
  {
    std::cout << "BIN MODE: Dynamic binning (bin fraction " << bin_fraction << " : " << int(ax.size()*bin_fraction) << " )" << std::endl;

    // bin construction
    vector<double> ax_sorted(ax), ay_sorted(ay), gz_sorted(gz);
    sort(ax_sorted.begin(), ax_sorted.end());
    sort(ay_sorted.begin(), ay_sorted.end());
    sort(gz_sorted.begin(), gz_sorted.end());

    int record_per_bin = int(ax.size()*bin_fraction);
    bin_number = int(1. / bin_fraction) + 1;

    std::cout << "AX range [ " << ax_min << " , " << ax_max << " ]\tbin_num: " << bin_number << std::endl;
    std::cout << "AY range [ " << ay_min << " , " << ay_max << " ]\tbin_num: " << bin_number << std::endl;
    std::cout << "GZ range [ " << gz_min << " , " << gz_max << " ]\tbin_num: " << bin_number << std::endl;

    vector<double> ax_pivot, ay_pivot, gz_pivot;
    for (int i = 1; i < bin_number; i++) {
      ax_pivot.push_back(ax_sorted[i*record_per_bin]);
      ay_pivot.push_back(ay_sorted[i*record_per_bin]);
      gz_pivot.push_back(gz_sorted[i*record_per_bin]);
    }
    ax_pivot.push_back(ax_max);
    ay_pivot.push_back(ay_max);
    gz_pivot.push_back(gz_max);

    //std::cout << "AX_PIVOT : "; for (auto p : ax_pivot) std::cout << p << "   "; std::cout << std::endl;
    //std::cout << "AY_PIVOT : "; for (auto p : ay_pivot) std::cout << p << "   "; std::cout << std::endl;
    //std::cout << "GZ_PIVOT : "; for (auto p : gz_pivot) std::cout << p << "   "; std::cout << std::endl;

    counter_ax_gz.resize(bin_number); for (auto &i : counter_ax_gz) i.resize(bin_number, 0);
    counter_ay_gz.resize(bin_number); for (auto &i : counter_ay_gz) i.resize(bin_number, 0);

    int ax_index, ay_index, gz_index;
    for (size_t i = 0; i < ax.size(); i++) {
      for (size_t j = 0; j < ax_pivot.size(); j++) if (ax[i] <= ax_pivot[j]) { ax_index = j; break; }
      for (size_t j = 0; j < ay_pivot.size(); j++) if (ay[i] <= ay_pivot[j]) { ay_index = j; break; }
      for (size_t j = 0; j < gz_pivot.size(); j++) if (gz[i] <= gz_pivot[j]) { gz_index = j; break; }

      //      std::cout << i << "\t" << ax_index << "\t" << ay_index << "\t" << gz_index << "\n";
      counter_ax_gz[ay_index][ax_index]++;
      counter_ay_gz[ay_index][gz_index]++;
    }

    //    int sum = 0;
    //    for (auto c : counter_ay_ax)
    //      for (auto e : c) sum += e;
    //    std::cout << "tot " << ax.size() << "\t" << sum << std::endl;

    break;
  }
  case FIXED_BIN:
  {
    std::cout << "BIN MODE: Fixed binning (bin number " << bin_number << " )" << std::endl;

    // bin_width calculation and output
    double ax_binw = (ax_max - ax_min) / double(bin_number - 1);
    double ay_binw = (ay_max - ay_min) / double(bin_number - 1);
    double gz_binw = (gz_max - gz_min) / double(bin_number - 1);
    std::cout << "AX range [ " << ax_min << " , " << ax_max << " ]\tbin_w: " << ax_binw << std::endl;
    std::cout << "AY range [ " << ay_min << " , " << ay_max << " ]\tbin_w: " << ay_binw << std::endl;
    std::cout << "GZ range [ " << gz_min << " , " << gz_max << " ]\tbin_w: " << gz_binw << std::endl;

    // counting frequencies
    counter_ax_gz.resize(bin_number); for (auto &i : counter_ax_gz) i.resize(bin_number, 0);   // improve this by maybe splitting the switch
    counter_ay_gz.resize(bin_number); for (auto &i : counter_ay_gz) i.resize(bin_number, 0);
    for (size_t i = 0; i < ay.size(); i++) {
      counter_ax_gz[size_t((ay[i] - ay_min) / ay_binw)][size_t((ax[i] - ax_min) / ax_binw)]++;
      counter_ay_gz[size_t((ay[i] - ay_min) / ay_binw)][size_t((gz[i] - gz_min) / gz_binw)]++;
    }

    break;
  }
  default:
    break;
  }

  // probability densities
  std::vector<double> p_ax(bin_number, 0), p_ay(bin_number, 0), p_gz(bin_number, 0);
  std::vector<std::vector<double>> p_ax_gz(bin_number, std::vector<double>(bin_number, 0.)), p_ay_gz(bin_number, std::vector<double>(bin_number, 0.));
  for (int i = 0; i < bin_number; i++) {
    for (int j = 0; j < bin_number; j++) {
      p_ax_gz[i][j] = counter_ax_gz[i][j] / double(ay.size());
      p_ay_gz[i][j] = counter_ay_gz[i][j] / double(ay.size());
      p_ax[i] += counter_ax_gz[j][i] / double(ay.size());
      p_ay[i] += counter_ay_gz[i][j] / double(ay.size());
      p_gz[i] += counter_ay_gz[j][i] / double(ay.size());
    }
  }

  //  double sum1 = 0, sum2 = 0, diff = 0;
  //  for (size_t i = 0; i < p_ay_a.size(); i++) {
  //    sum1 += p_ay_a[i];
  //    sum2 += p_ay_g[i];
  //    diff += abs(p_ay_a[i] - p_ay_g[i]);
  //  }
  //  std::cout << sum1 << "  " << sum2 << "  " << diff << std::endl;


    /*
    string outfile = input_file.substr(0, input_file.size() - 4) + "_count.txt";
    ofstream output(outfile);
    for (size_t i = 0; i < counter_ay_gz.size(); i++) {
      for (size_t j = 0; j < counter_ay_gz[i].size(); j++) {
        output << counter_ay_gz[i][j] << "\t";
      }
      output << "\t|\t" << counter_ay[i] << endl;
    }
    output << endl;
    for (size_t i = 0; i < counter_gz.size(); i++) {
      output << counter_gz[i] << "\t";
    }
    output.close();

    string densityfile = input_file.substr(0, input_file.size() - 4) + "_density.txt";
    output.open(densityfile);
    double I = 0;
    vector<double> p_ay(GRID_STEP, 0), p_gz(GRID_STEP, 0);
    vector<vector<double>> p_ay_gz(GRID_STEP, vector<double>(GRID_STEP, 0));
    for (size_t i = 0; i < counter_ay_gz.size(); i++) {
      p_ay[i] = counter_ay[i] / double(ay.size());
      for (size_t j = 0; j < counter_ay_gz[i].size(); j++) {
        p_ay_gz[i][j] = counter_ay_gz[i][j] / double(ay.size());
        p_gz[j] = counter_gz[j] / double(ay.size());

        output << i*ay_binw + ay_min << "\t" << j*gz_binw + gz_min << "\t" << p_ay_gz[i][j] << endl;
        I += ((p_ay_gz[i][j] == 0) ? 0 : p_ay_gz[i][j] * log(p_ay_gz[i][j] / (p_ay[i] * p_gz[j])));
        //cout << i << " " << j << " - " << p_ay_gz[i][j] << "(" << (p_ay_gz[i][j]==0) << ") " << p_ay[i] << " " << p_gz[j] << " " << log(p_ay_gz[i][j] / (p_ay[i] * p_gz[j])) << endl;
      }
    }
    output.close();

    string freqfile = input_file.substr(0, input_file.size() - 4) + "_freq.txt";
    output.open(freqfile);
    for (size_t i = 0; i < p_ay_gz.size(); i++) {
      for (size_t j = 0; j < p_ay_gz[i].size(); j++) {
        output << p_ay_gz[i][j] << "\t";
      }
      output << "\t|\t" << p_ay[i] << endl;
    }
    output << endl;
    for (size_t i = 0; i < counter_gz.size(); i++) {
      output << p_gz[i] << "\t";
    }
    output.close();
    */

  double me_ax_gz = mutual_entropy(p_ax_gz);
  double me_ay_gz = mutual_entropy(p_ay_gz);
  double e_ax = entropy(p_ax);
  double e_ay = entropy(p_ay);
  double e_gz = entropy(p_gz);

  std::cout << "AX entropy     = " << e_ax << endl;
  std::cout << "AY entropy     = " << e_ay << endl;
  std::cout << "GZ entropy     = " << e_gz << endl;
  std::cout << "AX-GZ mutual_e = " << me_ax_gz << endl;
  std::cout << "AY-GZ mutual_e = " << me_ay_gz << endl;
  std::cout << "Data samples   = " << ax.size() << endl;

  return 0;
}
