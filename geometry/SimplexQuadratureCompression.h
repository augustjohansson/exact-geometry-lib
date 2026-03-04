#pragma once

#include <cstddef>
#include <vector>
#include <utility>

namespace simpex
{

/// Compress a quadrature rule using the algorithm from:
///
///   A. Sommariva, M. Vianello,
///   "Compression of multivariate discrete measures and applications",
///   Numerical Functional Analysis and Optimization, 36(9), 2015.
///
/// \param[in,out] qr  Quadrature rule (flattened points and weights) to compress.
/// \param[in]     gdim Geometric dimension of the quadrature points.
/// \param[in]     quadrature_order Polynomial degree integrated exactly.
/// \returns The indices of retained quadrature points, or an empty vector if
///          no compression was performed.
std::vector<std::size_t>
compress_quadrature_rule(std::pair<std::vector<double>, std::vector<double>>& qr,
                         std::size_t gdim,
                         std::size_t quadrature_order);

} // namespace simpex
