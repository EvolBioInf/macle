#pragma once
#include <cstddef>
#include <cstdint>

/* a factorization of a string (LempelZiv or MatchLength) */
class Fact {
public:
  // Fact(Fact &other); //TODO: copy constructor, just to be sure? NRVO stuff
  ~Fact();
  void print() const;

  size_t *fact;     /* positions of factors */
  int64_t *prevOcc; /* LZ: prevOcc array (previous occurence of factor s_j) */
  size_t *lpf;      /* lpf in case of Lempel-Ziv, otherwise unset */
  size_t n;         /* number of factors */

  char const *str; /* string */
  size_t strLen;   /* string length */

  double cObs, cMax, cMin, cNor;
};

// length of factor
inline size_t factLen(Fact &f, size_t i) {
  if (i == 0)
    return f.fact[1];
  if (i == f.n - 1)
    return f.strLen - f.fact[f.n - 1];
  return f.fact[i + 1] - f.fact[i];
}
