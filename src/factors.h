#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include "config.h"

/* a factorization of a string (LempelZiv or MatchLength) */
struct Fact {
  void print() const;

  uint_vec fact;     /* positions of factors */
  uint_vec lpf;      /* lpf in case of Lempel-Ziv, otherwise empty */
  std::vector<int64_t> prevOcc; /* LZ: prevOcc array (previous occurence of factor s_j) */

  char const *str; /* string */
  size_t strLen;   /* string length */

  double cObs, cMax, cMin, cNor;
};

// length of factor
inline size_t factLen(Fact const &f, size_t i) {
  if (i == 0)
    return f.fact[1];
  if (i == f.fact.size() - 1)
    return f.strLen - f.fact[f.fact.size() - 1];
  return f.fact[i + 1] - f.fact[i];
}
