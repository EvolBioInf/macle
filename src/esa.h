/***** esa.h **************************************
 * Description: Header file for computation
 *   of enhance suffix array implemented
 *   in esa.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:17:08 2013
 **************************************************/
#pragma once

#include "prelude.h"
#include <divsufsort.h>

/* define data container */
typedef struct esa {
  saidx_t *sa;  /* suffix array */
  saidx_t *isa; /* inverse suffix array */
  int64_t *lcp; /* longest common prefix array */
  char const *str;    /* pointer to underlying string */
  size_t n;     /* length of sa and lcp */
} Esa;

Esa *getEsa(char const *seq, size_t n);
void freeEsa(Esa *esa);

void printEsa(Esa *esa);

int64_t *precomputeLcp(Esa *esa);
int64_t getLcp(Esa *esa, int64_t *lcptab, size_t sai, size_t saj);
