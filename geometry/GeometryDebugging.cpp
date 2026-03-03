// Copyright (C) 2016 Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2016-05-05
// Last changed: 2017-03-01

#include <sstream>
#include <iostream>

#include "GeometryDebugging.h"

using namespace dolfin;

// Plotting not initialized
bool GeometryDebugging::_initialized = false;

//-----------------------------------------------------------------------------
void GeometryDebugging::print(const Point& point)
{
  std::cout << "Point: " << point << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::print(const std::vector<Point>& simplex)
{
  std::cout << "Simplex:";
  for (const Point p : simplex)
    std::cout << " " << p;
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::print(const std::vector<Point>& simplex_0,
                              const std::vector<Point>& simplex_1)
{
  std::cout << "Simplex 0:";
  for (const Point p : simplex_0)
    std::cout << "-" << point2string(p);
  std::cout << std::endl;

  std::cout << "Simplex 1:";
  for (const Point p : simplex_1)
    std::cout << "-" << point2string(p);
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::plot(const Point& point)
{
  init_plot();

  std::cout << "# Plot point" << std::endl;
  std::cout << "ax.plot(" << simplex2string({point}) << ", 'x')" << std::endl;
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::plot(const std::vector<Point>& simplex)
{
  init_plot();

  std::cout << "# Plot simplex" << std::endl;
  if (simplex.size() >= 3)
  std::cout << "ax.plot_trisurf(" << simplex2string(simplex) << ")" << std::endl;
  else
    std::cout << "ax.plot(" << simplex2string(simplex) << ", marker='x')" << std::endl;
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::plot(const std::vector<Point>& simplex_0,
                             const std::vector<Point>& simplex_1)
{
  init_plot();

  std::cout << "# Plot simplex intersection" << std::endl;
  if (simplex_0.size() >= 3)
  std::cout << "ax.plot_trisurf(" << simplex2string(simplex_0) << ", color='r')" << std::endl;
  else
    std::cout << "ax.plot(" << simplex2string(simplex_0) << ", marker='o', color='r')" << std::endl;
  if (simplex_1.size() >= 3)
  std::cout << "ax.plot_trisurf(" << simplex2string(simplex_1) << ", color='b')" << std::endl;
  else
    std::cout << "ax.plot(" << simplex2string(simplex_1) << ", marker='o', color='b')" << std::endl;
  std::cout << std::endl;
}
//-----------------------------------------------------------------------------
void GeometryDebugging::init_plot()
{
  if (_initialized)
    return;

  std::cout << "# Initialize matplotlib 3D plotting" << std::endl;
  std::cout << "from mpl_toolkits.mplot3d import Axes3D" << std::endl;
  std::cout << "import matplotlib.pyplot as pl" << std::endl;
  std::cout << "ax = pl.figure().gca(projection='3d')" << std::endl;
  std::cout << "pl.ion(); pl.show()" << std::endl;
  std::cout << std::endl;
  std::cout << "# Note 1: Rotate/interact with figure to update plot." << std::endl;
  std::cout << "# Note 2: Use pl.cla() to clear figure between plots." << std::endl;
  std::cout << std::endl;

  _initialized = true;
}
//-----------------------------------------------------------------------------
std::string GeometryDebugging::point2string(const Point& p)
{
  std::stringstream s;
  s << "(" << p.x() << "," << p.y() << "," << p.z() << ")";
  return s.str();
}
//-----------------------------------------------------------------------------
std::string GeometryDebugging::simplex2string(const std::vector<Point>& simplex)
{
  std::size_t n = simplex.size();
  if (n == 0) return "";
  std::stringstream s;
  s << "[";
  for (std::size_t i = 0; i < n - 1; i++)
    s << simplex[i].x() << ",";
  s << simplex[n - 1].x() << "]";
  s << ",";
  s << "[";
  for (std::size_t i = 0; i < n - 1; i++)
    s << simplex[i].y() << ",";
  s << simplex[n - 1].y() << "]";
  s << ",";
  s << "[";
  for (std::size_t i = 0; i < n - 1; i++)
    s << simplex[i].z() << ",";
  s << simplex[n - 1].z() << "]";

  return s.str();
}
//-----------------------------------------------------------------------------
