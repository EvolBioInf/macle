/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#pragma once
#include <vector>
#include <divsufsort64.h>
#include "rmq.h"

/* define data container */
class Esa {
public:
  Esa(char const *seq, size_t n);
  RMQ precomputeLcp() const;
  int64_t getLcp(const RMQ &rmq, size_t sai, size_t saj) const;
  void print() const;

  std::vector<saidx64_t> sa;  /* suffix array */
  std::vector<saidx64_t> isa; /* inverse suffix array */
  std::vector<int64_t> lcp;   /* longest common prefix array */
  char const *str;            /* pointer to underlying string */
  size_t n;                   /* length of sa and lcp */
};
