// Copyright (C) 2016-2017 Anders Logg, August Johansson and Benjamin Kehlet
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

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <tuple>
#include <set>
#include "predicates.h"
#include "GeometryPredicates.h"
#include "GeometryTools.h"
#include "CollisionPredicates.h"
#include "IntersectionConstruction.h"
#include "ConvexTriangulation.h"
#include "CGALExactArithmetic.h"

using namespace simpex;

//-----------------------------------------------------------------------------
namespace
{
  // Create a unique list of points in the sense that |p-q| > tol in each dimension
  std::vector<Point>
  unique_points(const std::vector<Point>& input_points,
		std::size_t gdim,
		double tol)
  {
    std::vector<Point> points;

    for (std::size_t i = 0; i < input_points.size(); ++i)
    {
      bool unique = true;
      for (std::size_t j = i+1; j < input_points.size(); ++j)
      {
	std::size_t cnt = 0;
	for (std::size_t d = 0; d < gdim; ++d)
	{
	  if (std::abs(input_points[i][d] - input_points[j][d]) > tol)
	  {
	    cnt++;
	  }
	}
	if (cnt == 0)
	{
	  unique = false;
	  break;
	}
      }

      if (unique)
	points.push_back(input_points[i]);
    }

    return points;
  }
  //------------------------------------------------------------------------------
  // Check if q lies between p0 and p1. p0, p1 and q are assumed to be colinear
  inline bool is_between(const Point& p0, const Point& p1, const Point& q)
  {
    const double sqnorm = (p1-p0).squared_norm();
    return (p0-q).squared_norm() < sqnorm && (p1-q).squared_norm() < sqnorm;
  }
  //------------------------------------------------------------------------------
  // Return the edge pairs forming the convex hull of a set of coplanar points.
  std::vector<std::pair<std::size_t, std::size_t>>
  compute_convex_hull_planar(const std::vector<Point>& points)
  {
    Point normal = GeometryTools::cross_product(points[0], points[1], points[2]);
    normal /= normal.norm();

    std::vector<std::pair<std::size_t, std::size_t>> edges;

    for (std::size_t i = 0; i < points.size(); i++)
    {
      for (std::size_t j = i+1; j < points.size(); j++)
      {
        const Point r = points[i]+normal;

        double edge_orientation = 0;
        {
          std::size_t a = 0;
          while (edge_orientation == 0)
          {
            if (a != i && a != j)
            {
              edge_orientation = orient3d(points[i],
                                          points[j],
                                          r,
                                          points[a]);
            }
            a++;
          }
        }

        bool on_convex_hull = true;
	std::vector<std::size_t> colinear;
        for (std::size_t p = 0; p < points.size(); p++)
        {
          if (p != i && p != j)
          {
            const double orientation = orient3d(points[i],
                                                points[j],
                                                r,
                                                points[p]);

	    if (orientation == 0)
	    {
	      colinear.push_back(p);
	    }

            if (edge_orientation * orientation < 0)
            {
              on_convex_hull = false;
            }
          }
        }

        if (on_convex_hull)
        {
	  if (!colinear.empty())
	  {
	    bool is_linear_convex_hull = true;
	    for (std::size_t q : colinear)
	    {
	      if (!is_between(points[i], points[j], points[q]))
	      {
		is_linear_convex_hull = false;
		break;
	      }
	    }

	    if (is_linear_convex_hull)
	    {
	      edges.push_back(std::make_pair(i,  j));
	    }
	  }
	  else
	  {
	    edges.push_back(std::make_pair(i, j));
	  }
        }
      }
    }

    return edges;
  }
}

//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate(const std::vector<Point>& p,
                                 std::size_t gdim,
                                 std::size_t tdim)
{
  if (p.empty())
    return std::vector<std::vector<Point>>();

  if (tdim == 1)
  {
    return triangulate_1d(p, gdim);
  }
  else if (tdim == 2 && gdim == 2)
  {
    return triangulate_graham_scan_2d(p);
  }
  else if (tdim == 3 && gdim == 3)
  {
    return triangulate_graham_scan_3d(p);
  }

  throw std::runtime_error("Triangulation of polyhedron of topological dimension %u and geometric dimension %u not implemented");

  return std::vector<std::vector<Point>>();
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_1d(const std::vector<Point>& p,
				     std::size_t gdim)
{
  // A convex polyhedron of topological dimension 1 can not have more
  // than two points. If more, they must be collinear (more or
  // less). This can happen due to tolerances in
  // IntersectionConstruction::intersection_segment_segment_2d.

  if (gdim != 2)
  {
    throw std::runtime_error("Function is only implemented for gdim = 2");
  }

  const std::vector<Point> unique_p = unique_points(p, gdim, 3.0e-16);

  if (unique_p.size() == 2)
  {
    std::vector<std::vector<Point>> t { unique_p };
    return t;
  }
  else if (unique_p.size() < 2)
  {
    return std::vector<std::vector<Point>>();
  }
  else
  {
    // Here unique_p.size() > 2: points must be approximately collinear
    bool collinear = true;
    for (std::size_t i = 2; i < unique_p.size(); ++i)
    {
      const double o = orient2d(unique_p[0], unique_p[1], unique_p[i]);
      if (std::abs(o) > 1e-14)
      {
	collinear = false;
	break;
      }
    }

    assert(collinear);

    // Return extremal points
    const Point v = unique_p[1] - unique_p[0];
    std::vector<std::pair<double, std::size_t>> order;
    order.emplace_back(0.0, 0);
    order.emplace_back(1.0, 1);
    for (std::size_t i = 2; i < unique_p.size(); ++i)
      order.emplace_back(v.dot(unique_p[i]-unique_p[0]), i);

    std::sort(order.begin(), order.end());

    return {{ unique_p[order.front().second],
	  unique_p[order.back().second] }};
  }

}

//------------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_2d(const std::vector<Point>& input_points)
{
  assert(GeometryPredicates::is_finite(input_points));

  const std::size_t tdim = 2;
  const std::size_t gdim = 2;
  std::vector<Point> points = unique_points(input_points, gdim, 3.0e-16);

  if (points.size() < 3)
    return std::vector<std::vector<Point>>();

  if (points.size() == 3)
  {
    std::vector<std::vector<Point>> triangulation;
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }

  // Use the bottommost point (lowest y, break ties by lowest x) as the pivot.
  // This is the canonical Graham scan starting point.
  std::size_t pivot = 0;
  for (std::size_t m = 1; m < points.size(); ++m)
  {
    if (points[m].y() < points[pivot].y() ||
        (points[m].y() == points[pivot].y() && points[m].x() < points[pivot].x()))
      pivot = m;
  }
  std::swap(points[0], points[pivot]);

  // Compute center for angle reference
  Point pointscenter = points[0];
  for (std::size_t m = 1; m < points.size(); ++m)
    pointscenter += points[m];
  pointscenter /= points.size();

  const Point ref = points[0] - pointscenter;

  // Calculate and store angles
  std::vector<std::pair<double, std::size_t>> order;
  for (std::size_t m = 1; m < points.size(); ++m)
  {
    const double A = orient2d(pointscenter, points[0], points[m]);
    const Point s = points[m] - pointscenter;
    double alpha = std::atan2(A, s.dot(ref));
    if (alpha < 0)
      alpha += 2.0*3.141592653589793238462;
    order.emplace_back(alpha, m);
  }

  std::sort(order.begin(), order.end());

  std::vector<std::vector<Point>> triangulation;

  for (std::size_t m = 0; m < order.size()-1; ++m)
  {
    const std::vector<Point> tri {{ points[0],
	  points[order[m].second],
	  points[order[m + 1].second] }};
    if (!GeometryPredicates::is_degenerate(tri, tdim, gdim))
      triangulation.push_back(tri);
  }

  return triangulation;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::_triangulate_graham_scan_3d(const std::vector<Point>& input_points)
{
  assert(GeometryPredicates::is_finite(input_points));

  const std::size_t tdim = 3;
  const std::size_t gdim = 3;
  std::vector<Point> points = unique_points(input_points, gdim, 3.0e-16);

  std::vector<std::vector<Point>> triangulation;

  if (points.size() < 4)
  {
    return triangulation;
  }
  else if (points.size() == 4)
  {
    if (!GeometryPredicates::is_degenerate(points, tdim, gdim))
      triangulation.push_back(points);
    return triangulation;
  }
  else
  {
    Point polyhedroncenter(0,0,0);
    for (const Point& p : points)
      polyhedroncenter += p;
    polyhedroncenter /= points.size();

    std::set<std::tuple<std::size_t, std::size_t, std::size_t> > checked;

    for (std::size_t i = 0; i < points.size(); ++i)
    {
      for (std::size_t j = i+1; j < points.size(); ++j)
      {
        for (std::size_t k = j+1; k < points.size(); ++k)
        {
	  if (checked.emplace(std::make_tuple(i, j, k)).second)
	  {
            // Skip collinear triples using the exact cross-product (no division).
            if ((points[j]-points[i]).cross(points[k]-points[i]).squared_norm() == 0.0)
              continue;

	    bool on_convex_hull = true;
	    std::vector<std::size_t> coplanar = { i, j, k };
	    double previous_orientation = 0.0;
	    bool first = true;

	    for (std::size_t m = 0; m < points.size(); ++m)
	    {
	      if (m != i && m != j && m != k)
	      {
		const double orientation = orient3d(points[i],
						    points[j],
						    points[k],
						    points[m]);
		if (orientation == 0)
		  coplanar.push_back(m);
                else
                {
                  if (first)
                  {
                    previous_orientation = orientation;
                    first = false;
                  }
                  else
                  {
                    if (previous_orientation * orientation < 0)
                    {
                      on_convex_hull = false;
                    }
                  }
		}
	      }
	    }

	    if (on_convex_hull)
	    {
	      if (coplanar.size() == 3)
	      {
		std::vector<Point> cand = { points[i],
					    points[j],
					    points[k],
					    polyhedroncenter };
#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING
                if (cgal_tet_is_degenerate(cand))
                  throw std::runtime_error("tet is degenerate");
#endif
		if (!GeometryPredicates::is_degenerate(cand, tdim, gdim))
		  triangulation.push_back(cand);
	      }
	      else // At least four coplanar points
	      {
		std::vector<Point> coplanar_points;
		for (std::size_t idx : coplanar)
		  coplanar_points.push_back(points[idx]);

		std::vector<std::pair<std::size_t, std::size_t>> coplanar_convex_hull =
		  compute_convex_hull_planar(coplanar_points);

		Point coplanar_center(0,0,0);
		for (const Point& p : coplanar_points)
		  coplanar_center += p;
		coplanar_center /= coplanar_points.size();

		for (const std::pair<std::size_t, std::size_t>& edge : coplanar_convex_hull)
		{
		  const std::vector<Point> cand {{ polyhedroncenter,
			coplanar_center,
			coplanar_points[edge.first],
			coplanar_points[edge.second] }};
		  if (!GeometryPredicates::is_degenerate(cand, tdim, gdim))
		  {
		    triangulation.push_back(cand);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING
		    if (cgal_tet_is_degenerate(triangulation.back()))
		    {
		      throw std::runtime_error("tet is degenerate");
		    }

		    if (cgal_triangulation_overlap(triangulation))
		    {
		      throw std::runtime_error("now triangulation overlaps");
		    }
#endif
		  }

		  std::sort(coplanar.begin(), coplanar.end());

		  for (int coplanar_i = 0; coplanar_i < (int)coplanar.size()-2; coplanar_i++)
		  {
		    for (int coplanar_j = coplanar_i+1; coplanar_j < (int)coplanar.size()-1; coplanar_j++)
		    {
		      for (std::size_t coplanar_k = coplanar_j+1; coplanar_k < coplanar.size(); coplanar_k++)
		      {
			checked.emplace(std::make_tuple(coplanar[coplanar_i], coplanar[coplanar_j], coplanar[coplanar_k]));
		      }
                    }
                  }
                }
	      } // end coplanar.size() > 3
	    } // end on_convexhull
	  }
	}
      }
    }

    return triangulation;
  }
}
//-----------------------------------------------------------------------------
std::vector<std::vector<Point>>
ConvexTriangulation::triangulate_graham_scan_3d(const std::vector<Point>& pm)
{
  std::vector<std::vector<Point>> triangulation =
    _triangulate_graham_scan_3d(pm);

#ifdef DOLFIN_ENABLE_GEOMETRY_DEBUGGING

  if (cgal_triangulation_has_degenerate(triangulation))
    throw std::runtime_error("triangulation contains degenerate tetrahedron");

  if (cgal_triangulation_overlap(triangulation))
  {
    throw std::runtime_error("tetrahedrons overlap");
  }

  double volume = .0;
  for (const std::vector<Point>& tet : triangulation)
  {
    assert(tet.size() == 4);
    const double tet_volume = std::abs(orient3d(tet[0], tet[1], tet[2], tet[3]))/6.0;
    volume += tet_volume;
  }

  const double reference_volume = cgal_polyhedron_volume(pm);

  if (std::abs(volume - reference_volume) > 1e-14)
    throw std::runtime_error("computed volume %f, but reference volume is %f (diff %e)");

#endif
  return triangulation;
}
//-----------------------------------------------------------------------------
bool ConvexTriangulation::selfintersects(const std::vector<std::vector<Point>>& p)
{
  for (std::size_t i = 0; i < p.size(); i++)
  {
    for (std::size_t j = i+1; j < p.size(); j++)
    {
      assert(p[i].size() == p[j].size());
      if (p[i].size() == 4)
      {
	if (CollisionPredicates::collides_tetrahedron_tetrahedron_3d(p[i][0],
								     p[i][1],
								     p[i][2],
								     p[i][3],
								     p[j][0],
								     p[j][1],
								     p[j][2],
								     p[j][3]))
	{
	  auto intersection =
	    IntersectionConstruction::intersection_tetrahedron_tetrahedron_3d(p[i][0],
									      p[i][1],
									      p[i][2],
									      p[i][3],
									      p[j][0],
									      p[j][1],
									      p[j][2],
									      p[j][3]);
	  if (intersection.size() > 3)
	  {
	    for (std::size_t k = 3; k < intersection.size(); k++)
	    {
	      if (orient3d(intersection[0],
			   intersection[1],
			   intersection[2],
			   intersection[k]) != 0)
	      {
		return true;
	      }
	    }
	  }
	}
      }
      else if (p[i].size() == 3)
      {
	throw std::runtime_error("Not implemented");
      }
    }
  }

  return false;
}
//-----------------------------------------------------------------------------
