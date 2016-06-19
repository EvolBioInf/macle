#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>

#include "config.h"
#include "esa.h"

/* a match-length factorization of a string */
struct Fact {
  void print() const;
  size_t factLen(size_t i) const;

  uint_vec fact;     /* positions of factors */

  char const *str; /* string */
  size_t strLen;   /* string length */
};


void computeMLFact(Fact &fact, Esa const &esa);
