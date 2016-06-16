/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#pragma once
#include "config.h"
#include <divsufsort64.h>

/* define data container */
class Esa {
public:
  Esa(char const *seq, size_t n);
  void print() const;

  uint_vec sa;    /* suffix array */
  uint_vec isa;   /* inverse suffix array */
  uint_vec lcp;   /* longest common prefix array */
  char const *str;            /* pointer to underlying string */
  size_t n;                   /* length of sa and lcp */
};

void reduceEsa(Esa &esa);
