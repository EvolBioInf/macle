/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include "prelude.h"
#include "matchlength.h"
#include "eprintf.h"
#include "sequenceData.h"
#include "shulen.h"

// Fenwick-Tree (log(n) access to any prefix sums, n*log(n) construction)
void inc_ft(uint64_t *tree, size_t n, size_t i, uint64_t delta) {
  for (; i < n; i |= i + 1)
    tree[i] += delta;
}

int64_t sum_ft(uint64_t *tree, int64_t ind) {
  uint64_t sum = 0;
  while (ind >= 0) { // ind must be signed for it to work!!
    sum += tree[ind];
    ind &= ind + 1;
    ind--;
  }
  return sum;
}

uint64_t getsum_ft(uint64_t *tree, int64_t left, int64_t right) {
  return sum_ft(tree, right) - sum_ft(tree, left - 1);
}
//--------

uint64_t factorsFromTo(Fact *f, int64_t l, int64_t r) { return getsum_ft(f->lpf, l, r); }

Fact *computeMLFact(Esa *esa) {
  Fact *mlf = (Fact *)emalloc(sizeof(Fact));
  mlf->lpf = 0; // no lpf array in this factorization
  mlf->str = esa->str;
  mlf->strLen = esa->n;

  /* construct and fill array of match lengths */
  uint64_t *ml = (uint64_t *)emalloc(esa->n * sizeof(uint64_t));
  for (size_t i = 0; i < esa->n; i++) {
    ml[esa->sa[i]] = MAX(1, MAX(esa->lcp[i], esa->lcp[i + 1]));
  }

  /* compute observed number of match factors, store their positions */
  mlf->fact = (size_t *)emalloc(esa->n * sizeof(size_t));
  mlf->n = 0;
  size_t i = 0;
  while (i < esa->n) {
    mlf->fact[mlf->n++] = i;
    i += ml[i];
  }
  free(ml);
  return mlf;
}

// calculate number of factor start positions up to some given prefix length
// -> allows to easily get number of factor start positions in any interval
uint64_t *computeFactPrefixSum(Fact *f) {
  size_t n = f->strLen;
  /* compute observed number of match factors for every prefix */
  uint64_t *ft = (uint64_t *)ecalloc(n, sizeof(uint64_t));
  for (size_t i = 0; i < f->n; i++) {
    inc_ft(ft, n, f->fact[i], 1);
  }
  return ft;
}

/* mlComplexity: calculate match length factors */
Fact *mlComplexity(Esa *esa, double gc) {
  size_t n = esa->n;
  /* construct and fill array of match lengths, calculate factors */
  Fact *mlf = computeMLFact(esa);
  /* compute observed number of match factors for every prefix */
  uint64_t *fl = computeFactPrefixSum(mlf);
  mlf->lpf = fl; // use lpf field to store prefix sums for ML factors

  // calculation of global complexity value (for reference):
  // TODO: refactor this out

  double c = (double)(sum_ft(fl, n - 1)); // total number of factors
  double l = n;
  mlf->cObs = (c - 1) / l; // TODO: why c-1?
  mlf->cMin = 2. / l;

  /* TODO: need this?
  // cludge to prevent infinite loop &  non-sensical result from expShulen
  if(gc == 1.0)
    gc = 0.9;
  */
  double esl = expShulen(gc, n);
  mlf->cMax = 1. / (esl - 1.);
  // scale value to relative amount between min. and estimated max.
  mlf->cNor = (mlf->cObs - mlf->cMin) / (mlf->cMax - mlf->cMin);
  //----

  return mlf;
}
