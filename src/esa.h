/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#pragma once

#include "prelude.h"

/* define data container */
typedef struct esa {
  uint64_t *sa; /* suffix array */
  int64_t *lcp; /* longest common prefix array */
  char *str;    /* pointer to underlying string */
  size_t n;     /* length of sa and lcp */
} Esa;

Esa *getEsa(char *seq, size_t n);
void freeEsa(Esa *esa);
