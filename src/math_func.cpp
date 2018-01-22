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

#include <cstddef>
#include <vector>

#include "libbbmutils/bbmutils.h"
#include "math_func.h"
#include "params.h"

std::vector<double> Integrate(std::vector<double>& x, std::vector<double>& f, double F0)
{
  std::vector<double> F;
  F.push_back(F0);
  double dx, df, sum = F[0];
  for (size_t i = 1; i < x.size(); i++) {
    dx = x[i] - x[i - 1];
    df = .5 * (f[i] + f[i - 1]);
    sum += dx * df;
    F.push_back(sum);
  }
  return F;
}

std::vector<std::vector<double>> rotate_inertial(std::vector<std::vector<double>> data, MAT3D rotation)
{
  std::vector<std::vector<double>> data_r;

  for (size_t i = 0; i < data.size(); i++) {
    std::vector<double> row;
    row.push_back(data[i][0]); // timestamp
    row.push_back(data[i][1]); // speed

    VEC6D ag;
    set_vec6d(&ag, data[i][AX_INDEX], data[i][AY_INDEX], data[i][AZ_INDEX],
        data[i][GX_INDEX], data[i][GY_INDEX], data[i][GZ_INDEX]);
    VEC6D ag_r;
    rotate_vec6d(&ag_r, &rotation, &ag);
    row.push_back(ag_r.a.x);
    row.push_back(ag_r.a.y);
    row.push_back(ag_r.a.z);
    row.push_back(ag_r.g.x);
    row.push_back(ag_r.g.y);
    row.push_back(ag_r.g.z);

    for (size_t j = 8; j < data[i].size(); j++) {
      row.push_back(data[i][j]);
    }
    data_r.push_back(row);
  }

  return data_r;
}

std::vector<double> forward_derivative(std::vector<double> F, std::vector<double> x)
{
  std::vector<double> f;
  for (size_t i = 0; i < F.size(); i++) {
    if (i == (F.size() - 1)) {
      f.push_back(f[i - 1]);
    } else {
      f.push_back((F[i + 1] - F[i]) / (x[i + 1] - x[i]));
    }
  }
  return f;
}
