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
    static inline std::size_t choose(std::size_t n, std::size_t k)
    {
      if (k > n) return 0;
      if (k == 0 || k == n) return 1;
      if (k > n - k) k = n - k;
      long double res = 1.0L;
      for (std::size_t i = 1; i <= k; ++i)
        res = res * static_cast<long double>(n - k + i) / static_cast<long double>(i);
      return static_cast<std::size_t>(std::llround(res));
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
  compress(std::pair<std::vector<double>, std::vector<double>>& qr,
	   std::size_t gdim,
	   std::size_t quadrature_order,
	   CompressionMethod method)
  {
    switch (method)
      {
      case CompressionMethod::legacy:
	return compress_legacy(qr, gdim, quadrature_order);

      case CompressionMethod::nnls:
	return compress_nnls(qr, gdim, quadrature_order);

      case CompressionMethod::irls:
	return compress_irls(qr, gdim, quadrature_order);
      }

    throw std::runtime_error("Unknown compression method");
  }

  // ---------------------------------------------------------------------------
  std::vector<std::size_t>
  compress_legacy(std::pair<std::vector<double>, std::vector<double>>& qr,
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

static inline Eigen::VectorXd nnls_lawson_hanson(
						 const Eigen::MatrixXd& A,
						 const Eigen::VectorXd& b,
						 double tol = 1e-12,
						 int max_outer = 5000,
						 int max_inner = 5000)
{
  const int m = static_cast<int>(A.rows());
  const int n = static_cast<int>(A.cols());

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd w = A.transpose() * (b - A * x);

  std::vector<bool> P(n, false);
  std::vector<bool> Z(n, true);

  auto max_w_in_Z = [&]() -> std::pair<int,double> {
    int t = -1;
    double maxw = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i)
      if (Z[i] && w[i] > maxw) { maxw = w[i]; t = i; }
    return {t, maxw};
  };

  int outer = 0;
  while (outer++ < max_outer)
    {
      auto [t, maxw] = max_w_in_Z();
      if (t < 0 || maxw <= tol) break;

      P[t] = true; Z[t] = false;

      int inner = 0;
      while (inner++ < max_inner)
        {
	  std::vector<int> pidx;
	  pidx.reserve(n);
	  for (int i = 0; i < n; ++i) if (P[i]) pidx.push_back(i);

	  const int p = static_cast<int>(pidx.size());
	  Eigen::MatrixXd A_P(m, p);
	  for (int j = 0; j < p; ++j) A_P.col(j) = A.col(pidx[j]);

	  Eigen::VectorXd z = A_P.colPivHouseholderQr().solve(b);

	  Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
	  for (int j = 0; j < p; ++j) s[pidx[j]] = z[j];

	  bool all_pos = true;
	  for (int j = 0; j < p; ++j)
	    if (s[pidx[j]] <= tol) { all_pos = false; break; }

	  if (all_pos)
            {
	      x = s;
	      break;
            }

	  double alpha = std::numeric_limits<double>::infinity();
	  for (int j = 0; j < p; ++j)
            {
	      int idx = pidx[j];
	      if (s[idx] <= tol)
                {
		  double denom = x[idx] - s[idx];
		  if (denom > 0.0)
		    alpha = std::min(alpha, x[idx] / denom);
                }
            }
	  if (!std::isfinite(alpha)) alpha = 0.0;

	  x = x + alpha * (s - x);

	  for (int i = 0; i < n; ++i)
            {
	      if (P[i] && x[i] <= tol)
                {
		  x[i] = 0.0;
		  P[i] = false;
		  Z[i] = true;
                }
            }
        }

      w = A.transpose() * (b - A * x);
    }

  for (int i = 0; i < x.size(); ++i) if (x[i] < 0.0) x[i] = 0.0;
  return x;
}

static inline Eigen::VectorXd irls_l1(
				      const Eigen::MatrixXd& A,
				      const Eigen::VectorXd& b,
				      double lambda,
				      int iters = 30,
				      double eps = 1e-12)
{
  const Eigen::MatrixXd AtA = A.transpose() * A;
  const Eigen::VectorXd Atb = A.transpose() * b;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
  for (int k = 0; k < iters; ++k)
    {
      Eigen::VectorXd w = (x.array().abs() + eps).inverse().matrix();
      Eigen::MatrixXd M = AtA + lambda * w.asDiagonal();
      x = M.ldlt().solve(Atb);
    }
  return x;
}

static inline std::vector<std::size_t> overwrite_with_indices(
							      std::pair<std::vector<double>, std::vector<double>>& qr,
							      const std::pair<std::vector<double>, std::vector<double>>& qr_input,
							      const Eigen::VectorXd& w_new,
							      std::size_t gdim,
							      std::size_t K_bound,
							      double tol_keep)
{
  const std::size_t M = qr_input.second.size();

  std::vector<std::size_t> idx;
  idx.reserve(M);
  for (Eigen::Index i = 0; i < w_new.size(); ++i)
    if (std::abs(w_new[i]) > tol_keep)
      idx.push_back(static_cast<std::size_t>(i));

  // If solver returns too many nonzeros, keep the largest |w|
  if (idx.size() > K_bound)
    {
      std::vector<std::pair<double,std::size_t>> mag;
      mag.reserve(idx.size());
      for (auto i : idx) mag.emplace_back(std::abs(w_new[static_cast<Eigen::Index>(i)]), i);

      std::nth_element(mag.begin(), mag.begin() + static_cast<long>(K_bound), mag.end(),
		       [](const auto& a, const auto& b){ return a.first > b.first; });
      mag.resize(K_bound);

      idx.clear();
      idx.reserve(K_bound);
      for (auto& p : mag) idx.push_back(p.second);
      std::sort(idx.begin(), idx.end());
    }

  qr.first.resize(gdim * idx.size());
  qr.second.resize(idx.size());

  for (std::size_t i = 0; i < idx.size(); ++i)
    {
      const std::size_t src = idx[i];
      for (std::size_t d = 0; d < gdim; ++d)
	qr.first[gdim * i + d] = qr_input.first[gdim * src + d];
      qr.second[i] = w_new[static_cast<Eigen::Index>(src)];
    }

  return idx;
}

static inline Eigen::VectorXd nnls_lawson_hanson(
						 const Eigen::MatrixXd& A,
						 const Eigen::VectorXd& b,
						 double tol = 1e-12,
						 int max_outer = 5000,
						 int max_inner = 5000)
{
  const int m = static_cast<int>(A.rows());
  const int n = static_cast<int>(A.cols());

  Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd w = A.transpose() * (b - A * x);

  std::vector<bool> P(n, false);
  std::vector<bool> Z(n, true);

  auto max_w_in_Z = [&]() -> std::pair<int,double> {
    int t = -1;
    double maxw = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i)
      if (Z[i] && w[i] > maxw) { maxw = w[i]; t = i; }
    return {t, maxw};
  };

  int outer = 0;
  while (outer++ < max_outer)
    {
      auto [t, maxw] = max_w_in_Z();
      if (t < 0 || maxw <= tol) break;

      P[t] = true; Z[t] = false;

      int inner = 0;
      while (inner++ < max_inner)
        {
	  std::vector<int> pidx;
	  pidx.reserve(n);
	  for (int i = 0; i < n; ++i) if (P[i]) pidx.push_back(i);

	  const int p = static_cast<int>(pidx.size());
	  Eigen::MatrixXd A_P(m, p);
	  for (int j = 0; j < p; ++j) A_P.col(j) = A.col(pidx[j]);

	  Eigen::VectorXd z = A_P.colPivHouseholderQr().solve(b);

	  Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
	  for (int j = 0; j < p; ++j) s[pidx[j]] = z[j];

	  bool all_pos = true;
	  for (int j = 0; j < p; ++j)
	    if (s[pidx[j]] <= tol) { all_pos = false; break; }

	  if (all_pos)
            {
	      x = s;
	      break;
            }

	  double alpha = std::numeric_limits<double>::infinity();
	  for (int j = 0; j < p; ++j)
            {
	      int idx = pidx[j];
	      if (s[idx] <= tol)
                {
		  double denom = x[idx] - s[idx];
		  if (denom > 0.0)
		    alpha = std::min(alpha, x[idx] / denom);
                }
            }
	  if (!std::isfinite(alpha)) alpha = 0.0;

	  x = x + alpha * (s - x);

	  for (int i = 0; i < n; ++i)
            {
	      if (P[i] && x[i] <= tol)
                {
		  x[i] = 0.0;
		  P[i] = false;
		  Z[i] = true;
                }
            }
        }

      w = A.transpose() * (b - A * x);
    }

  for (int i = 0; i < x.size(); ++i) if (x[i] < 0.0) x[i] = 0.0;
  return x;
}

static inline Eigen::VectorXd irls_l1(
				      const Eigen::MatrixXd& A,
				      const Eigen::VectorXd& b,
				      double lambda,
				      int iters = 30,
				      double eps = 1e-12)
{
  const Eigen::MatrixXd AtA = A.transpose() * A;
  const Eigen::VectorXd Atb = A.transpose() * b;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(A.cols());
  for (int k = 0; k < iters; ++k)
    {
      Eigen::VectorXd w = (x.array().abs() + eps).inverse().matrix();
      Eigen::MatrixXd M = AtA + lambda * w.asDiagonal();
      x = M.ldlt().solve(Atb);
    }
  return x;
}


// Positive-weight compression (TETRAFREEQ-like)
static inline std::vector<std::size_t>
compress_nnls(std::pair<std::vector<double>, std::vector<double>>& qr,
              std::size_t gdim,
              std::size_t quadrature_order)
{
  const std::size_t N = quadrature_order;
  const std::size_t K = choose(N + gdim, gdim);
  const std::size_t M = qr.second.size();

  if (M <= K) return {};

  const auto qr_input = qr;
  const Eigen::MatrixXd V = chebyshev_vandermonde_matrix(qr_input, gdim, N);

  if (static_cast<std::size_t>(V.rows()) != K || static_cast<std::size_t>(V.cols()) != M)
    throw std::runtime_error("Chebyshev_Vandermonde_matrix returned wrong shape.");

  Eigen::Map<const Eigen::VectorXd> w_base(qr_input.second.data(),
					   static_cast<Eigen::Index>(M));
  const Eigen::VectorXd b = V * w_base;

  // NNLS: min ||V w - b||, w >= 0
  Eigen::VectorXd w_new = nnls_lawson_hanson(V, b, /*tol=*/1e-12);

  const double maxabs = w_new.cwiseAbs().maxCoeff();
  const double tol_keep = std::max(1e-15, 1e-12 * maxabs);

  return overwrite_with_indices(qr, qr_input, w_new, gdim, K, tol_keep);
}

// Signed-weight sparse compression (CheapQ-compatible)
static inline std::vector<std::size_t>
compress_signed_irls(std::pair<std::vector<double>, std::vector<double>>& qr,
                     std::size_t gdim,
                     std::size_t quadrature_order)
{
  const std::size_t N = quadrature_order;
  const std::size_t K = choose(N + gdim, gdim);
  const std::size_t M = qr.second.size();

  if (M <= K) return {};

  const auto qr_input = qr;
  const Eigen::MatrixXd V = chebyshev_vandermonde_matrix(qr_input, gdim, N);

  if (static_cast<std::size_t>(V.rows()) != K || static_cast<std::size_t>(V.cols()) != M)
    throw std::runtime_error("Chebyshev_Vandermonde_matrix returned wrong shape.");

  Eigen::Map<const Eigen::VectorXd> w_base(qr_input.second.data(),
					   static_cast<Eigen::Index>(M));
  const Eigen::VectorXd b = V * w_base;

  // IRLS parameters: for N=2 or 4, K is small (2D: 6/15; 3D: 10/35)
  // so we can push toward exactness while encouraging sparsity.
  const double lambda = 1e-12; // tune if needed
  Eigen::VectorXd w_new = irls_l1(V, b, lambda, /*iters=*/40, /*eps=*/1e-12);

  const double maxabs = w_new.cwiseAbs().maxCoeff();
  const double tol_keep = std::max(1e-15, 1e-12 * maxabs);

  return overwrite_with_indices(qr, qr_input, w_new, gdim, K, tol_keep);
}
