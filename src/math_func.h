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

#include <vector>

#include "math_struct.h"

std::vector<double> Integrate(std::vector<double> &x, std::vector<double> &f, double F0 = 0);

std::vector< std::vector<double> > rotate_inertial(std::vector< std::vector<double> > data, MAT3D rotation);

std::vector<double> forward_derivative(std::vector<double> F, std::vector<double> x);