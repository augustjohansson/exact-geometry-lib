#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#include <Eigen/Dense>

#include "SimplexQuadratureCompression.h"

namespace simpex
{

namespace
{

// ---------------------------------------------------------------------------
// Helper: number of combinations n choose k
// ---------------------------------------------------------------------------
static std::size_t choose(std::size_t n, std::size_t k)
{
  if (k == 0) return 1;
  return (n * choose(n - 1, k - 1)) / k;
}

// ---------------------------------------------------------------------------
// Helper: Chebyshev polynomials of the first kind T_0 ... T_N evaluated at x
// ---------------------------------------------------------------------------
static std::vector<Eigen::VectorXd>
chebyshev_polynomial(const Eigen::VectorXd& x, std::size_t N)
{
  // T_0(x) = 1
  // T_1(x) = x
  // T_k(x) = 2 x T_{k-1}(x) - T_{k-2}(x)

  std::vector<Eigen::VectorXd> T(N + 1, Eigen::VectorXd(x.size()));
  T[0] = Eigen::VectorXd::Ones(x.size());

  if (N >= 1)
  {
    T[1] = x;
    for (std::size_t k = 2; k <= N; ++k)
      T[k] = 2.0 * x.cwiseProduct(T[k - 1]) - T[k - 2];
  }

  return T;
}

// ---------------------------------------------------------------------------
// Helper: graded lexicographic ordering of multi-indices up to total degree N
// ---------------------------------------------------------------------------
static std::vector<std::vector<std::size_t>>
grlex(std::size_t gdim, std::size_t N)
{
  const std::size_t n_rows = choose(N + gdim, gdim);
  std::vector<std::vector<std::size_t>> P(n_rows, std::vector<std::size_t>(gdim));

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
    assert(false && "grlex only implemented for gdim = 2 or 3");
  }

  return P;
}

// ---------------------------------------------------------------------------
// Helper: Chebyshev–Vandermonde matrix for the given quadrature rule
// ---------------------------------------------------------------------------
static Eigen::MatrixXd
chebyshev_vandermonde_matrix(
    const std::pair<std::vector<double>, std::vector<double>>& qr,
    std::size_t gdim,
    std::size_t N)
{
  const auto npts = static_cast<Eigen::Index>(qr.second.size());
  std::vector<std::vector<Eigen::VectorXd>> T(gdim);
  Eigen::VectorXd x(npts);

  for (std::size_t d = 0; d < gdim; ++d)
  {
    for (std::size_t i = 0; i < static_cast<std::size_t>(npts); ++i)
      x(static_cast<Eigen::Index>(i)) = qr.first[i * gdim + d];

    const double xmin = x.minCoeff();
    const double xmax = x.maxCoeff();
    const double hx   = xmax - xmin;
    const Eigen::VectorXd xmap =
        (2.0 / hx) * x - ((xmin + xmax) / hx) * Eigen::VectorXd::Ones(x.size());

    T[d] = chebyshev_polynomial(xmap, N);
  }

  const std::vector<std::vector<std::size_t>> P = grlex(gdim, N);
  const std::size_t n_cols = P.size();

  Eigen::MatrixXd V = Eigen::MatrixXd::Ones(npts, static_cast<Eigen::Index>(n_cols));

  for (std::size_t i = 1; i < n_cols; ++i)
  {
    Eigen::VectorXd V_col = T[0][P[i][0]];
    for (std::size_t d2 = 1; d2 < gdim; ++d2)
      V_col = V_col.cwiseProduct(T[d2][P[i][d2]]);
    V.col(static_cast<Eigen::Index>(i)) = V_col;
  }

  return V;
}

} // anonymous namespace

// ---------------------------------------------------------------------------

std::vector<std::size_t>
compress_quadrature_rule(std::pair<std::vector<double>, std::vector<double>>& qr,
                         std::size_t gdim,
                         std::size_t quadrature_order)
{
  const std::size_t N = quadrature_order;
  const std::size_t N_min = choose(N + gdim, gdim);

  if (qr.second.size() <= N_min)
    return {};

  const auto qr_input = qr;

  const Eigen::MatrixXd V = chebyshev_vandermonde_matrix(qr_input, gdim, N);

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(V);
  Eigen::MatrixXd Q = QR.householderQ();
  Q *= Eigen::MatrixXd::Identity(V.rows(), std::min(V.rows(), V.cols()));
  Q.transposeInPlace();

  Eigen::Map<const Eigen::VectorXd> w_base(
      qr_input.second.data(),
      static_cast<Eigen::Index>(qr_input.second.size()));

  const Eigen::VectorXd nu    = Q * w_base;
  const Eigen::VectorXd w_new = Q.colPivHouseholderQr().solve(nu);

  std::vector<std::size_t> indices;
  for (Eigen::Index i = 0; i < w_new.size(); ++i)
    if (std::abs(w_new[i]) > 0.0)
      indices.push_back(static_cast<std::size_t>(i));

  assert(indices.size() <= N_min);

  qr.first.resize(gdim * indices.size());
  qr.second.resize(indices.size());

  for (std::size_t i = 0; i < indices.size(); ++i)
  {
    for (std::size_t d = 0; d < gdim; ++d)
      qr.first[gdim * i + d] = qr_input.first[gdim * indices[i] + d];
    qr.second[i] = w_new[static_cast<Eigen::Index>(indices[i])];
  }

  return indices;
}

} // namespace simpex
