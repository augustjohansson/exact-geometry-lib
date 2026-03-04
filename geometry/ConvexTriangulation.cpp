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


struct ProjPt
{
  double u, v;
  std::size_t idx;
};

static inline int dominant_drop_axis_from_normal(const Point& n) noexcept
{
  const double ax = std::abs(n.x());
  const double ay = std::abs(n.y());
  const double az = std::abs(n.z());
  if (ax >= ay && ax >= az) return 0; // drop x -> use (y,z)
  if (ay >= ax && ay >= az) return 1; // drop y -> use (x,z)
  return 2;                           // drop z -> use (x,y)
}

static inline void project_drop_axis(const Point& p, const Point& origin,
                                     int drop, double& u, double& v) noexcept
{
  const double x = p.x() - origin.x();
  const double y = p.y() - origin.y();
  const double z = p.z() - origin.z();
  if (drop == 0)      { u = y; v = z; }
  else if (drop == 1) { u = x; v = z; }
  else                { u = x; v = y; }
}

// Find a “good” normal direction by selecting the max-area triangle.
// Returns false if all points are collinear/identical (no plane).
static bool find_stable_normal(const std::vector<Point>& P, Point& N_out,
                               std::size_t& anchor_out)
{
  const std::size_t n = P.size();
  double best = 0.0;
  Point bestN;
  std::size_t bestA = 0;

  for (std::size_t a = 0; a < n; ++a)
    for (std::size_t b = a + 1; b < n; ++b)
      for (std::size_t c = b + 1; c < n; ++c)
      {
        Point N = GeometryTools::cross_product(P[a], P[b], P[c]); // (b-a)x(c-a)
        const double s = N.squared_norm();
        if (s > best)
        {
          best = s;
          bestN = N;
          bestA = a;
        }
      }

  if (best == 0.0)
    return false;

  N_out = bestN;     // NOTE: not normalized
  anchor_out = bestA;
  return true;
}

// Hybrid turn: orient2d fast; if 0, use orient3d(a,b,a+N,c) fallback.
static inline double turn_sign_hybrid(const ProjPt& A, const ProjPt& B, const ProjPt& C,
                                      const std::vector<Point>& P,
                                      const Point& N)
{
  const double a2[2] = {A.u, A.v};
  const double b2[2] = {B.u, B.v};
  const double c2[2] = {C.u, C.v};
  const double o2 = _orient2d(a2, b2, c2);
  if (o2 != 0.0) return o2;

  // Fallback: robust in-plane orientation using 3D predicate
  const Point& a3 = P[A.idx];
  const Point& b3 = P[B.idx];
  const Point& c3 = P[C.idx];
  const Point r = a3 + N;
  return orient3d(a3, b3, r, c3);
}

static std::vector<std::size_t>
convex_hull_planar_indices_robust(const std::vector<Point>& points)
{
  const std::size_t n = points.size();
  if (n == 0) return {};
  if (n == 1) return {0};

  Point N;
  std::size_t anchor = 0;
  if (!find_stable_normal(points, N, anchor))
  {
    // 1D hull: choose endpoints along axis with largest spread
    double minx = points[0].x(), maxx = points[0].x();
    double miny = points[0].y(), maxy = points[0].y();
    double minz = points[0].z(), maxz = points[0].z();
    for (std::size_t i = 1; i < n; ++i)
    {
      minx = std::min(minx, points[i].x()); maxx = std::max(maxx, points[i].x());
      miny = std::min(miny, points[i].y()); maxy = std::max(maxy, points[i].y());
      minz = std::min(minz, points[i].z()); maxz = std::max(maxz, points[i].z());
    }
    const double sx = maxx - minx, sy = maxy - miny, sz = maxz - minz;
    int axis = (sx >= sy && sx >= sz) ? 0 : (sy >= sx && sy >= sz) ? 1 : 2;

    auto coord = [&](const Point& p) -> double {
      return (axis == 0) ? p.x() : (axis == 1) ? p.y() : p.z();
    };

    std::size_t lo = 0, hi = 0;
    for (std::size_t i = 1; i < n; ++i)
    {
      if (coord(points[i]) < coord(points[lo])) lo = i;
      if (coord(points[i]) > coord(points[hi])) hi = i;
    }
    if (lo == hi) return {lo};
    return {lo, hi};
  }

  const int drop = dominant_drop_axis_from_normal(N);
  const Point origin = points[anchor];

  // Project
  std::vector<ProjPt> P2;
  P2.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
  {
    double u, v;
    project_drop_axis(points[i], origin, drop, u, v);
    P2.push_back({u, v, i});
  }

  // Sort by (u,v)
  std::sort(P2.begin(), P2.end(), [](const ProjPt& a, const ProjPt& b) {
    if (a.u < b.u) return true;
    if (a.u > b.u) return false;
    if (a.v < b.v) return true;
    if (a.v > b.v) return false;
    return a.idx < b.idx;
  });

  // Remove exact duplicates in projection
  P2.erase(std::unique(P2.begin(), P2.end(), [](const ProjPt& a, const ProjPt& b) {
    return a.u == b.u && a.v == b.v;
  }), P2.end());

  if (P2.size() == 1) return {P2[0].idx};
  if (P2.size() == 2) return {P2[0].idx, P2[1].idx};

  // Andrew monotone chain, pop on <=0 to discard collinear boundary points
  std::vector<ProjPt> lower, upper;
  lower.reserve(P2.size());
  upper.reserve(P2.size());

  for (const auto& p : P2)
  {
    while (lower.size() >= 2)
    {
      const auto& A = lower[lower.size() - 2];
      const auto& B = lower[lower.size() - 1];
      const double o = turn_sign_hybrid(A, B, p, points, N);
      if (o <= 0.0) lower.pop_back();
      else break;
    }
    lower.push_back(p);
  }

  for (std::size_t k = P2.size(); k-- > 0; )
  {
    const auto& p = P2[k];
    while (upper.size() >= 2)
    {
      const auto& A = upper[upper.size() - 2];
      const auto& B = upper[upper.size() - 1];
      const double o = turn_sign_hybrid(A, B, p, points, N);
      if (o <= 0.0) upper.pop_back();
      else break;
    }
    upper.push_back(p);
  }

  lower.pop_back();
  upper.pop_back();

  std::vector<std::size_t> hull;
  hull.reserve(lower.size() + upper.size());
  for (const auto& p : lower) hull.push_back(p.idx);
  for (const auto& p : upper) hull.push_back(p.idx);

  return hull;
}

std::vector<std::pair<std::size_t, std::size_t>>
compute_convex_hull_planar(const std::vector<Point>& points)
{
  std::vector<std::pair<std::size_t, std::size_t>> edges;

  const auto hull = convex_hull_planar_indices_robust(points);
  if (hull.size() < 2) return edges;

  edges.reserve(hull.size());
  for (std::size_t i = 0; i < hull.size(); ++i)
  {
    const std::size_t a = hull[i];
    const std::size_t b = hull[(i + 1) % hull.size()];
    edges.emplace_back(a, b);
  }

  return edges;
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
