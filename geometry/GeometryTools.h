// Copyright (C) 2017 Anders Logg
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
// First added:  2017-02-11
// Last changed: 2017-02-17

#ifndef __GEOMETRY_TOOLS_H
#define __GEOMETRY_TOOLS_H

#include <cassert>
#include "predicates.h"
#include "Point.h"

namespace simpex
{

  /// This class provides useful tools (functions) for computational geometry.

  class GeometryTools
  {
  public:

    /// Compute numerically stable cross product (a - c) x (b - c)
    static inline Point cross_product(const Point& a, const Point& b, const Point& c)
    {
      // See Shewchuk Lecture Notes on Geometric Robustness
      double ayz[2] = {a.y(), a.z()};
      double byz[2] = {b.y(), b.z()};
      double cyz[2] = {c.y(), c.z()};
      double azx[2] = {a.z(), a.x()};
      double bzx[2] = {b.z(), b.x()};
      double czx[2] = {c.z(), c.x()};
      double axy[2] = {a.x(), a.y()};
      double bxy[2] = {b.x(), b.y()};
      double cxy[2] = {c.x(), c.y()};
      return Point (_orient2d(ayz, byz, cyz),
                    _orient2d(azx, bzx, czx),
                    _orient2d(axy, bxy, cxy));
    }

    /// Compute determinant of 3 x 3 matrix defined by vectors, ab, dc, ec
    inline double determinant(const Point& ab, const Point& dc, const Point& ec)
    {
      // const double a = ab.x(), b = ab.y(), c = ab.z();
      // const double d = dc.x(), e = dc.y(), f = dc.z();
      // const double g = ec.x(), h = ec.y(), i = ec.z();
      // return a * (e * i - f * h)
      //      + b * (f * g - d * i)
      //      + c * (d * h - e * g);

      const double fast = ab.x() * (dc.y() * ec.z() - dc.z() * ec.y())
	+ ab.y() * (dc.z() * ec.x() - dc.x() * ec.z())
	+ ab.z() * (dc.x() * ec.y() - dc.y() * ec.x());

      // If far from zero, return fast.
      // Threshold can be tuned; this is a conservative relative-ish guard.
      const double s = std::abs(fast);
      if (s > 1e-18) return fast;

      static const Point O(0,0,0);
      return orient3d(ab, dc, ec, O);
    }

    /// Compute major (largest) axis of vector (2D)
    static inline std::size_t major_axis_2d(const Point& v)
    {
      return (std::abs(v.x()) >= std::abs(v.y()) ? 0 : 1);
    }

    /// Compute major (largest) axis of vector (3D)
    static inline std::size_t major_axis_3d(const Point& v)
    {
      const double vx = std::abs(v.x());
      const double vy = std::abs(v.y());
      const double vz = std::abs(v.z());
      if (vx >= vy && vx >= vz)
        return 0;
      if (vy >= vz)
        return 1;
      return 2;
    }

    /// Project point to axis (2D)
    static inline double project_to_axis_2d(const Point& p, std::size_t axis)
    {
      assert(axis <= 1);
      return p[axis];
    }

    /// Project point to plane (3D)
    static inline Point project_to_plane_3d(const Point& p, std::size_t axis)
    {
      assert(axis <= 2);
      switch (axis)
	{
	case 0: return Point(p.y(), p.z());
	case 1: return Point(p.x(), p.z());
	case 2: return Point(p.x(), p.y());
	}
      return p;
    }

    /// Check whether x in [a, b]
    static inline bool contains(double a, double b, double x)
    {
      return a <= x and x <= b;
    }

    /// Check whether x in (a, b)
    static inline bool contains_strict(double a, double b, double x)
    {
      return a < x and x < b;
    }

    /// Compute segment length (1D, 2D, or 3D)
    static inline double segment_length(const Point& a, const Point& b, std::size_t gdim)
    {
      const Point e = b - a;
      if (gdim == 1) return std::abs(e.x());
      if (gdim == 2) return std::sqrt(e.x()*e.x() + e.y()*e.y());
      return e.norm();
    }

    /// Compute triangle area (2D in-plane or 3D embedded)
    static inline double triangle_area(const Point& a, const Point& b, const Point& c,
                                        std::size_t gdim)
    {
      if (gdim == 2) return 0.5 * std::abs(orient2d(a, b, c));
      const Point u = b - a;
      const Point v = c - a;
      return 0.5 * u.cross(v).norm();
    }

    /// Compute tetrahedron volume
    static inline double tetra_volume(const Point& a, const Point& b,
                                       const Point& c, const Point& d)
    {
      return std::abs(orient3d(a, b, c, d)) / 6.0;
    }

  };

}

#endif
