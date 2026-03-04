// Copyright (C) 2014 August Johansson and Anders Logg
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

#include <cassert>
#include <cmath>
#include <stdexcept>
#include "GeometryTools.h"
#include "SimplexQuadrature.h"
#include "SimplexQuadratureTables.h"
#include "predicates.h"

using namespace simpex;

// ---------------------------------------------------------------------------
// File-scope helpers
// ---------------------------------------------------------------------------
static void write_point(double* out, std::size_t gdim, const Point& p)
{
  out[0] = p.x();
  if (gdim > 1) out[1] = p.y();
  if (gdim > 2) out[2] = p.z();
}

//-----------------------------------------------------------------------------
SimplexQuadrature::SimplexQuadrature(std::size_t tdim, std::size_t order)
{
  // Create and store quadrature rule for reference simplex
  switch (tdim)
    {
    case 1:
      setup_qr_reference_interval(order);
      break;
    case 2:
      setup_qr_reference_triangle(order);
      break;
    case 3:
      setup_qr_reference_tetrahedron(order);
      break;
    default:
      throw std::runtime_error("Only implemented for topological dimension 1, 2, 3");
    }

}
//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule(const std::vector<Point>& coordinates,
					   std::size_t gdim) const
{
  std::size_t tdim = coordinates.size() - 1;

  switch (tdim)
    {
    case 0:
      // FIXME: should we return empty qr or should we have detected this earlier?
      break;
    case 1:
      return compute_quadrature_rule_interval(coordinates, gdim);
      break;
    case 2:
      return compute_quadrature_rule_triangle(coordinates, gdim);
      break;
    case 3:
      return compute_quadrature_rule_tetrahedron(coordinates, gdim);
      break;
    default:
      throw std::runtime_error("Only implemented for topological dimension 1, 2, 3");
    };

  std::pair<std::vector<double>, std::vector<double>> quadrature_rule;
  return quadrature_rule;
}
//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_interval(const std::vector<Point>& X,
                                                    std::size_t gdim) const
{
  if (X.size() != 2) throw std::runtime_error("interval: need 2 vertices");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * gdim);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e  = X[1] - X[0];

  const double L = GeometryTools::segment_length(X[0], X[1], gdim);
  const double ref = 1.0;                              // |[0,1]|
  const double scale = L / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const Point p = x0 + e * r;
      write_point(&pts[i*gdim], gdim, p);
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_triangle(const std::vector<Point>& X,
                                                    std::size_t gdim) const
{
  if (X.size() != 3) throw std::runtime_error("triangle: need 3 vertices");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * gdim);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e1 = X[1] - X[0];
  const Point  e2 = X[2] - X[0];

  const double A = GeometryTools::triangle_area(X[0], X[1], X[2], gdim);
  const double ref = 0.5; // reference triangle area
  const double scale = A / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const double s = _p[i][1];
      const Point p = x0 + e1 * r + e2 * s;
      write_point(&pts[i*gdim], gdim, p);
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
std::pair<std::vector<double>, std::vector<double>>
SimplexQuadrature::compute_quadrature_rule_tetrahedron(const std::vector<Point>& X,
                                                       std::size_t gdim) const
{
  if (X.size() != 4) throw std::runtime_error("tetrahedron: need 4 vertices");
  if (gdim != 3) throw std::runtime_error("tetrahedron: requires gdim=3");
  const std::size_t n = _p.size();

  std::vector<double> pts(n * 3);
  std::vector<double> w(n);

  const Point& x0 = X[0];
  const Point  e1 = X[1] - X[0];
  const Point  e2 = X[2] - X[0];
  const Point  e3 = X[3] - X[0];

  const double V = GeometryTools::tetra_volume(X[0], X[1], X[2], X[3]);
  const double ref = 1.0 / 6.0;
  const double scale = V / ref;

  for (std::size_t i = 0; i < n; ++i)
    {
      const double r = _p[i][0];
      const double s = _p[i][1];
      const double t = _p[i][2];
      const Point p = x0 + e1 * r + e2 * s + e3 * t;
      pts[i*3 + 0] = p.x();
      pts[i*3 + 1] = p.y();
      pts[i*3 + 2] = p.z();
      w[i] = _w[i] * scale;
    }

  return {std::move(pts), std::move(w)};
}

//-----------------------------------------------------------------------------
std::vector<std::size_t>
SimplexQuadrature::compress(std::pair<std::vector<double>, std::vector<double>>& qr,
			    std::size_t gdim,
			    std::size_t quadrature_order)
{
  // Polynomial degree N that can be integrated exactly using the
  // qr_base
  const std::size_t N = quadrature_order;

  // By construction the compressed quadrature rule will not have more
  // than choose(N + gdim, gdim) points.
  const std::size_t N_compressed_min = choose(N + gdim, gdim);
  if (qr.second.size() <= N_compressed_min)
    {
      // We cannot improve this rule. Return empty vector
      return std::vector<std::size_t>();
    }

  // Copy the input qr since we'll overwrite the input
  const std::pair<std::vector<double>, std::vector<double>> qr_input = qr;

  // Create Vandermonde-type matrix using a basis of Chebyshev
  // polynomials of the first kind
  const Eigen::MatrixXd V = Chebyshev_Vandermonde_matrix(qr_input, gdim, N);

  // A QR decomposition selects the subset of N columns (geometrically
  // the N columns with same volume as spanned by all M columns).
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(V);
  Eigen::MatrixXd Q = QR.householderQ();

  // We do not need the full Q matrix but only what's known as the
  // "economy size" decomposition
  Q *= Eigen::MatrixXd::Identity(V.rows(), std::min(V.rows(), V.cols()));

  // We'll use Q^T
  Q.transposeInPlace();

  // Compute weights in the new basis
  Eigen::Map<const Eigen::VectorXd> w_base(qr_input.second.data(),
					   static_cast<Eigen::Index>(qr_input.second.size()));
  const Eigen::VectorXd nu = Q*w_base;

  // Compute new weights
  const Eigen::VectorXd w_new = Q.colPivHouseholderQr().solve(nu);

  // Construct new qr using the non-zero weights. First find the
  // indices for these weights.
  std::vector<std::size_t> indices;
  for (Eigen::Index i = 0; i < w_new.size(); ++i)
    {
      if (std::abs(w_new[i]) > 0.0)
	indices.push_back(static_cast<std::size_t>(i));
    }

  // Resize qr and overwrite the points and weights
  assert(indices.size() <= N_compressed_min);
  qr.first.resize(gdim*indices.size());
  qr.second.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); ++i)
    {
      // Save points
      for (std::size_t d = 0; d < gdim; ++d)
	qr.first[gdim*i + d] = qr_input.first[gdim*indices[i] + d];

      // Save weights
      qr.second[i] = w_new[static_cast<Eigen::Index>(indices[i])];
    }

  // Return indices. These are useful for mapping additional data, for
  // example the normals.
  return indices;
}
//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_interval(std::size_t order)
{
  const auto rr = tables::interval_rule_legendre(order);
  _p.assign(rr.n, std::vector<double>(1));
  _w.assign(rr.w, rr.w + rr.n);
  for (std::size_t i = 0; i < rr.n; ++i)
    _p[i][0] = rr.pts[3*i + 0]; // r
}

//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_triangle(std::size_t order)
{
  if (order > 20) throw std::runtime_error("Triangle quadrature supported up to order 20");
  const auto rr = tables::triangle_rule_dunavant(order);
  _p.assign(rr.n, std::vector<double>(2));
  _w.assign(rr.w, rr.w + rr.n);
  for (std::size_t i = 0; i < rr.n; ++i)
    {
      _p[i][0] = rr.pts[3*i + 0]; // r
      _p[i][1] = rr.pts[3*i + 1]; // s
    }
}

//-----------------------------------------------------------------------------
void SimplexQuadrature::setup_qr_reference_tetrahedron(std::size_t order)
{
  if (order > 8) throw std::runtime_error("Tetrahedron quadrature supported up to order 8");
  const auto rr = tables::tetra_rule_keast(order);
  _p.assign(rr.n, std::vector<double>(3));
  _w.assign(rr.w, rr.w + rr.n);
  for (std::size_t i = 0; i < rr.n; ++i)
    {
      _p[i][0] = rr.pts[3*i + 0]; // r
      _p[i][1] = rr.pts[3*i + 1]; // s
      _p[i][2] = rr.pts[3*i + 2]; // t
    }
}
//-----------------------------------------------------------------------------
Eigen::MatrixXd SimplexQuadrature::Chebyshev_Vandermonde_matrix
(const std::pair<std::vector<double>, std::vector<double>>& qr,
 std::size_t gdim,
 std::size_t N)
{
  // Create the Chebyshev basis matrix for each dimension separately
  std::vector<std::vector<Eigen::VectorXd>> T(gdim);
  Eigen::VectorXd x(static_cast<Eigen::Index>(qr.second.size()));

  for (std::size_t d = 0; d < gdim; ++d)
    {
      // Extract coordinates in one dimension
      for (std::size_t i = 0; i < qr.second.size(); ++i)
	x(static_cast<Eigen::Index>(i)) = qr.first[i*gdim + d];

      // Map points to [-1, 1]
      const double xmin = x.minCoeff();
      const double xmax = x.maxCoeff();
      const double hx = xmax - xmin;
      const Eigen::VectorXd xmap = (2./hx) * x - ((xmin+xmax)/hx) * Eigen::VectorXd::Ones(x.size());

      // Evaluate the basis
      T[d] = Chebyshev_polynomial(xmap, N);
    }

  // Find the order of the polynomials in graded lexicographic
  // ordering
  const std::vector<std::vector<std::size_t>> P = grlex(gdim, N);

  // Setup the Vandermonde type matrix
  const std::size_t n_cols = P.size();
  Eigen::MatrixXd V(Eigen::MatrixXd::Ones(static_cast<Eigen::Index>(qr.second.size()),
                                          static_cast<Eigen::Index>(n_cols)));

  // The first column is always [1, 1, ..., 1], hence we can start from 1
  for (std::size_t i = 1; i < n_cols; ++i)
    {
      // Pick out the correct order of polynomial from the P
      // matrix. Start with dimension 0.
      Eigen::VectorXd V_col = T[0][P[i][0]];

      // Do a .* style multiplication for each other dimension
      for (std::size_t d2 = 1; d2 < gdim; ++d2)
	{
	  V_col = V_col.cwiseProduct(T[d2][P[i][d2]]);
	}
      V.col(static_cast<Eigen::Index>(i)) = V_col;
    }

  return V;
}
//-----------------------------------------------------------------------------
std::vector<Eigen::VectorXd>
SimplexQuadrature::Chebyshev_polynomial(const Eigen::VectorXd& x,
					std::size_t N)
{
  // Create Chebyshev polynomial of the first kind of order N on [-1, 1]
  //
  // T_0(x) = 1
  // T_1(x) = x
  // T_{k}(x) = 2 * x * T_{k-1}(x) - T_{k-2}(x), k = 2, ..., N

  // Store in a matrix T such that (T)_ij = T_i(x_j). We don't use
  // Eigen::MatrixXd, because we want to slice out the rows later.
  std::vector<Eigen::VectorXd> T(N + 1, Eigen::VectorXd(x.size()));

  // Reccurence construction
  T[0] = Eigen::VectorXd::Ones(x.size());

  if (N >= 1)
    {
      T[1] = x;
      for (std::size_t k = 2; k <= N; ++k)
	T[k] = 2*x.cwiseProduct(T[k-1]) - T[k-2];
    }

  return T;
}
//-----------------------------------------------------------------------------
std::vector<std::vector<std::size_t>>
SimplexQuadrature::grlex(std::size_t gdim, std::size_t N)
{
  // Generate a matrix with numbers in graded lexicographic ordering,
  // i.e. if N = 3 and dim = 2, P should be
  // P = [ 0 0
  //       0 1
  //       1 0
  //       0 2
  //       1 1
  //       2 0
  //       0 3
  //       1 2
  //       2 1
  //       3 0 ]

  const std::size_t n_rows = choose(N + gdim, gdim);
  std::vector<std::vector<std::size_t>> P(n_rows, std::vector<std::size_t>(gdim));

  // FIXME: Make this a dimension independent loop
  switch (gdim)
    {
    case 2:
      for (std::size_t sum = 0, row = 0; sum <= N; ++sum)
	for (std::size_t xi = 0; xi <= N; ++xi)
	  for (std::size_t yi = 0; yi <= N; ++yi)
	    if (xi + yi == sum)
	      {
		P[row][0] = xi;
		P[row][1] = yi;
		row++;
	      }
      break;
    case 3:
      for (std::size_t sum = 0, row = 0; sum <= N; ++sum)
	for (std::size_t xi = 0; xi <= N; ++xi)
	  for (std::size_t yi = 0; yi <= N; ++yi)
	    for (std::size_t zi = 0; zi <= N; ++zi)
	      if (xi + yi + zi == sum)
		{
		  P[row][0] = xi;
		  P[row][1] = yi;
		  P[row][2] = zi;
		  row++;
		}
      break;
    default:
      assert(false);
    }

  return P;
}
//-----------------------------------------------------------------------------
std::size_t
SimplexQuadrature::choose(std::size_t n,
			  std::size_t k)
{
  // Compute the number of combinations n over k
  if (k == 0)
    return 1;
  return (n * choose(n - 1, k - 1)) / k;
}
//-----------------------------------------------------------------------------
