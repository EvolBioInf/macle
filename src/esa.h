/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#pragma once
#include <divsufsort64.h>

#include <sdsl/int_vector.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/rmq_support.hpp>

/* define data container */
class Esa {
public:
  Esa(char const *seq, size_t n);
  sdsl::rmq_succinct_sct<> precomputeLcp() const;
  int64_t getLcp(sdsl::rmq_succinct_sct<> const &rmq, size_t sai, size_t saj) const;
  void print() const;

  sdsl::int_vector<VECBIT> sa;    /* suffix array */
  sdsl::int_vector<VECBIT> isa;   /* inverse suffix array */
  sdsl::int_vector<VECBIT> lcp;   /* longest common prefix array */
  char const *str;            /* pointer to underlying string */
  size_t n;                   /* length of sa and lcp */
};

void reduceEsa(Esa &esa);
