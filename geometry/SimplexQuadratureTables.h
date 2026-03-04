#pragma once
#include <cstddef>

namespace simpex::tables {

  struct RuleRef
  {
    std::size_t n;
    const double* pts; // packed as (r,s,t) triples; unused entries are 0
    const double* w;
  };

  // interval rules are generated on [0,1] as (r,0,0), sum(w)=1
  RuleRef interval_rule_legendre(std::size_t order);

  // triangle rules on reference triangle (0,0),(1,0),(0,1), sum(w)=0.5
  RuleRef triangle_rule_dunavant(std::size_t order); // order<=20

  // tetra rules on reference tetra (0,0,0),(1,0,0),(0,1,0),(0,0,1), sum(w)=1/6
  RuleRef tetra_rule_keast(std::size_t order); // order<=8

} // namespace simpex::tables
