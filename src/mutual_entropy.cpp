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


#define MAJOR_VERSION       3
#define MINOR_VERSION       5

enum bin_mode { DYNAMIC_BIN, FIXED_BIN };

#include "params.h"

#include <iostream>
#include <iomanip>
#include <complex>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "io_lib.hpp"

using namespace std;

//////// Math functions
double entropy(const vector<double>& p) {
  double e = 0.;
  for (auto& pi : p) if (pi > numeric_limits<double>::epsilon()) e += pi*log(pi);
  return -e;
}


double mutual_entropy(const vector<vector<double>>& pxy) {
  vector<double> px(pxy.size(), 0.), py(pxy[0].size(), 0.);
  for (size_t i = 0; i < pxy.size(); i++) {
    for (size_t j = 0; j < pxy[i].size(); j++) {
      px[i] += pxy[i][j];
      py[j] += pxy[i][j];
    }
  }

  double me = 0.;
  for (size_t i = 0; i < pxy.size(); i++) {
    for (size_t j = 0; j < pxy[i].size(); j++) {
      if (!(pxy[i][j] > numeric_limits<double>::epsilon())) continue;
      me += pxy[i][j] * log(pxy[i][j] / (px[i] * py[j]));
    }
  }
  return me;
}


void usage(char * progname) {
  std::vector<std::string> tokens;
  boost::split(tokens, progname, boost::is_any_of("/\\"));
  cerr << "Usage: " << tokens.back() << " [BIN_MODE] [SHIFT_MODE] path/to/file" << endl;
  cerr << "       [BIN_MODE] (mandatory)" << endl;
  cerr << "                 -bf <bin_fraction> : DYNAMIC binning, represents the fraction of data in each bin (e.g. 0.05 corresponds to 5%)." << endl;
  cerr << "                 -bn <bin_number>   : FIXED binning, represents the number of bins for each set of data." << endl;
  cerr << "       [SHIFT_MODE] (optional)" << endl;
  cerr << "                 -s <shift>             : SINGLE shift, represents the shift in indices between GZ and AX,AY (negative values allowed, safe)" << endl;
  cerr << "                 -ss <min> <max> <incr> : SCAN shift, safe and non UNIX sed compliant" << endl;
  cerr << "       path/to/data/file must be a PHYSYCOM INERTIAL STANDARD" << endl << endl;
}


int main(int argc, char **argv) {
  cout << "Mutual Entropy Calculator v" << MAJOR_VERSION << "." << MINOR_VERSION << endl;
  string input_file;
  double bin_fraction;
  int bin_number;
  vector<int> index_shifts;
  char mode;
  if (argc > 3) {
    for (int i = 0; i < argc; i++) {
      if (string(argv[i]) == "-bf") {
        mode = DYNAMIC_BIN;
        try { bin_fraction = stod(string(argv[++i])); }
        catch (exception &e) {
          cerr << "EXCEPTION: " << e.what() << endl;
          usage(argv[0]);
          exit(-4);
        }
      }
      if (string(argv[i]) == "-bn") {
        mode = FIXED_BIN;
        try { bin_number = stoi(string(argv[++i])); }
        catch (exception &e) {
          cerr << "EXCEPTION: " << e.what() << endl;
          usage(argv[0]);
          exit(-5);
        }
      }
      if (string(argv[i]) == "-s") {
        try { index_shifts.push_back(stoi(string(argv[++i]))); }
        catch (exception &e) {
          cerr << "EXCEPTION: " << e.what() << endl;
          usage(argv[0]);
          exit(-6);
        }
      }
      if (string(argv[i]) == "-ss") {
        try {
          int start_shift = stoi(string(argv[++i]));
          int last_shift = stoi(string(argv[++i]));
          int scan_shift = stoi(string(argv[++i]));
          if (last_shift < start_shift) {
            int temp = start_shift;
            start_shift = last_shift;
            last_shift = temp;
          }
          if (scan_shift < 1) {
            cerr << "Invalid scan shift, reset to 1" << endl;
            scan_shift = 1;
          }
          for (int k = start_shift; k <= last_shift; k += scan_shift) index_shifts.push_back(k);
        }
        catch (exception &e) {
          cerr << "EXCEPTION: " << e.what() << endl;
          usage(argv[0]);
          exit(-7);
        }
      }
    }
    input_file = argv[argc - 1];
    if (input_file.substr(0, 2) == ".\\") input_file = input_file.substr(2, input_file.size() - 2);
  }
  else {
    cerr << "ERROR: Wrong command line parameters. Read usage and relaunch properly." << endl;
    usage(argv[0]);
    exit(-3);
  }


  // data parsing, conversion, storage
  vector< vector<string> > file_tokens = Read_from_file(input_file);
  vector< vector<double> > data = tokens_to_double(file_tokens);

  if (!data.size()) {
    cerr << "ERROR: Empty file " << input_file << endl;
    exit(-4);
  }

  if (file_tokens.size() != data.size()) {
    cout << "WARNING: original record " << file_tokens.size() << ", used " << data.size() << " ( " << int(100.0*(double)data.size()/file_tokens.size()) << " % )" << endl;
    cout << "WARNING: skipped " << file_tokens.size() - data.size() << " lines due to NaN presence." << endl;
  }

  vector<double> me_ax_gz(index_shifts.size(), std::numeric_limits<double>::signaling_NaN());
  vector<double> me_ay_gz(index_shifts.size(), std::numeric_limits<double>::signaling_NaN());
  vector<double> e_ax(index_shifts.size(), std::numeric_limits<double>::signaling_NaN());
  vector<double> e_ay(index_shifts.size(), std::numeric_limits<double>::signaling_NaN());
  vector<double> e_gz(index_shifts.size(), std::numeric_limits<double>::signaling_NaN());
  vector<int> bin_pop(index_shifts.size(), std::numeric_limits<int>::signaling_NaN());
  
  vector<double> ax_, ay_, gz_;
  for (auto line : data) {
    ax_.push_back(line[AX_INDEX]);
    ay_.push_back(line[AY_INDEX]);
    gz_.push_back(line[GZ_INDEX]);
  }

#pragma omp parallel for
  for (int n = 0; n < index_shifts.size(); n++) {
    if (abs(index_shifts[n]) > (int)ax_.size() / 2) {
      //cout << "Index shift " << index_shifts[n] << " is too big, the upper bound is: " << ax_.size() / 2 << endl;
      continue;
    }

    vector<double> ax;
    vector<double> ay;
    vector<double> gz;
    for (size_t i = 0; i < ax_.size() - abs(index_shifts[n]); i++) {
      if (index_shifts[n] >= 0) {
        ax.push_back(ax_[i]);
        ay.push_back(ay_[i]);
        gz.push_back(gz_[i + index_shifts[n]]);
      }
      else {
        ax.push_back(ax_[i - index_shifts[n]]);
        ay.push_back(ay_[i - index_shifts[n]]);
        gz.push_back(gz_[i]);
      }
    }

    // ranges
    auto   ax_mm = minmax_element(ax.begin(), ax.end());
    double ax_min = *ax_mm.first;
    double ax_max = *ax_mm.second;

    auto   ay_mm = minmax_element(ay.begin(), ay.end());
    double ay_min = *ay_mm.first;
    double ay_max = *ay_mm.second;

    auto   gz_mm = minmax_element(gz.begin(), gz.end());
    double gz_min = *gz_mm.first;
    double gz_max = *gz_mm.second;

    // binning
    vector<vector<int>> counter_ax_gz, counter_ay_gz;
    int record_per_bin = 1;
    switch (mode) {
    case DYNAMIC_BIN:
    {
      // bin construction
      vector<double> ax_sorted(ax), ay_sorted(ay), gz_sorted(gz);
      sort(ax_sorted.begin(), ax_sorted.end());
      sort(ay_sorted.begin(), ay_sorted.end());
      sort(gz_sorted.begin(), gz_sorted.end());

      record_per_bin = int(ax.size()*bin_fraction);
      bin_number = int(1. / bin_fraction) + 1;

      vector<double> ax_pivot, ay_pivot, gz_pivot;
      for (int i = 1; i < bin_number; i++) {
        ax_pivot.push_back(ax_sorted[i*record_per_bin]);
        ay_pivot.push_back(ay_sorted[i*record_per_bin]);
        gz_pivot.push_back(gz_sorted[i*record_per_bin]);
      }
      ax_pivot.push_back(ax_max);
      ay_pivot.push_back(ay_max);
      gz_pivot.push_back(gz_max);
      counter_ax_gz.resize(bin_number); for (auto &i : counter_ax_gz) i.resize(bin_number, 0);
      counter_ay_gz.resize(bin_number); for (auto &i : counter_ay_gz) i.resize(bin_number, 0);

      size_t ax_index, ay_index, gz_index;
      for (size_t i = 0; i < ax.size(); i++) {
        for (size_t j = 0; j < ax_pivot.size(); j++) if (ax[i] <= ax_pivot[j]) { ax_index = j; break; }
        for (size_t j = 0; j < ay_pivot.size(); j++) if (ay[i] <= ay_pivot[j]) { ay_index = j; break; }
        for (size_t j = 0; j < gz_pivot.size(); j++) if (gz[i] <= gz_pivot[j]) { gz_index = j; break; }

        counter_ax_gz[ay_index][ax_index]++;
        counter_ay_gz[ay_index][gz_index]++;
      }
      break;
    }
    case FIXED_BIN:
    {
      record_per_bin = int(ax.size()) / bin_number;

      // bin_width calculation
      double ax_binw = (ax_max - ax_min) / double(bin_number - 1);
      double ay_binw = (ay_max - ay_min) / double(bin_number - 1);
      double gz_binw = (gz_max - gz_min) / double(bin_number - 1);

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
    vector<double> p_ax(bin_number, 0.0);
    vector<double> p_ay(bin_number, 0.0);
    vector<double> p_gz(bin_number, 0.0);
    vector<vector<double>> p_ax_gz(bin_number, vector<double>(bin_number, 0.));
    vector<vector<double>> p_ay_gz(bin_number, vector<double>(bin_number, 0.));
    double renorm = 1.0 / double(ay.size());
    for (int i = 0; i < bin_number; i++) {
      for (int j = 0; j < bin_number; j++) {
        p_ax_gz[i][j] = counter_ax_gz[i][j] * renorm;
        p_ay_gz[i][j] = counter_ay_gz[i][j] * renorm;
        p_ax[i] += counter_ax_gz[j][i] * renorm;
        p_ay[i] += counter_ay_gz[i][j] * renorm;
        p_gz[i] += counter_ay_gz[j][i] * renorm;
      }
    }

    me_ax_gz[n] = mutual_entropy(p_ax_gz);
    me_ay_gz[n] = mutual_entropy(p_ay_gz);
    e_ax[n] = entropy(p_ax);
    e_ay[n] = entropy(p_ay);
    e_gz[n] = entropy(p_gz);
    bin_pop[n] = record_per_bin;

  } // omp parallel for ends here

  // dumping results
  std::string results_filename = input_file.substr(0, input_file.size() - 4) + "_entropy.txt";
  std::string gnuplot_filename = input_file.substr(0, input_file.size() - 4) + "_entropy.plt";
  std::string plot_filename = input_file.substr(0, input_file.size() - 4) + "_entropy.png";
  std::string escaped_filename = boost::replace_all_copy(results_filename, "_", "\\_");

  ofstream results(results_filename);
  results << "## Shift scan @ " << ((mode == DYNAMIC_BIN) ? "bin_fraction : " + to_string(bin_fraction) : "bin_number : " + to_string(bin_number)) << " # Data samples : " << ax_.size() << endl;
  results << "# shift # Entropy AX # Entropy AY # Entropy GZ # Mutual AX-GZ # Mutual AY-GZ # bin population #" << endl;
  for (int i = 0; i < index_shifts.size(); i++) {
    if (bin_pop[i] == bin_pop[i]) results << index_shifts[i] << "\t" << e_ax[i] << "\t" << e_ay[i] << "\t" << e_gz[i] << "\t" << me_ax_gz[i] << "\t" << me_ay_gz[i] << "\t" << bin_pop[i] << endl;
    else continue; // skip results still initialized to NaN (meaning they were not calculated in the algoritm)
  }
  results.close();

  ofstream gnuplot(gnuplot_filename);
  gnuplot << R"(#!/gnuplot
FILE_IN=')" << results_filename << R"('
FILE_OUT=')" << plot_filename << R"('
set terminal pngcairo dashed size 1280,720 enhanced font 'Verdana,10'
set output FILE_OUT
# Styles
linew = 1.2
set style line  21 lc rgb '#0072bd' lt 7 lw linew  # blue
set style line  22 lc rgb '#d95319' lt 7 lw linew  # orange
set style line  23 lc rgb '#77ac30' lt 7 lw linew  # green
set style line  24 lc rgb '#a2142f' lt 7 lw linew  # red
set style line 102 lc rgb '#d6d7d9' lt 1 lw 1      # gray
# Grid
set grid xtics ytics back ls 102
# Titles
set title 'Mutual Entropy: )" << escaped_filename << R"('
set xlabel 'Index Shift'
set ylabel 'Entropy'
set y2label 'Bin Population'
set ytics nomirror
set y2tics
# Plot
plot FILE_IN u 1:5 w lines ls 21 t 'ME(a_x , g_z)' axes x1y1,\
     FILE_IN u 1:6 w lines ls 24 t 'ME(a_y , g_z)' axes x1y1,\
     FILE_IN u 1:7 w lines ls 23 t 'Bin Pop'       axes x1y2
)";
  gnuplot.close();
  return 0;
}


