/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include "prelude.h"
#include "matchlength.h"
#include "eprintf.h"
#include "shulen.h"

uint64_t factorsFromTo(Fact *f, int64_t l, int64_t r) { return f->lpf[r]-(l?f->lpf[l-1]:0); }

Fact *computeMLFact(Esa *esa) {
  Fact *mlf = (Fact *)emalloc(sizeof(Fact));
  mlf->lpf = NULL;     // no lpf array in this factorization
  mlf->prevOcc = NULL; // no prevOcc too
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
void computeFactPrefixSum(Fact *f) {
  size_t n = f->strLen;
  /* compute observed number of match factors for every prefix */
  uint64_t *ps = (uint64_t *)emalloc(n * sizeof(uint64_t));
  size_t nextfact = 1;
  ps[0]=1;
  for (size_t i = 1; i < n; i++) {
    ps[i] = ps[i-1];
    if (nextfact<f->n && i==f->fact[nextfact]) {
      ps[i]++;
      nextfact++;
    }
  }
  f->lpf = ps; // use lpf field to store prefix sums for ML factors
}

/* mlComplexity: calculate match length factors */
Fact *mlComplexity(Esa *esa, double gc) {
  size_t n = esa->n;
  /* construct and fill array of match lengths, calculate factors */
  Fact *mlf = computeMLFact(esa);
  /* compute observed number of match factors for every prefix */
  computeFactPrefixSum(mlf);

  // calculation of global complexity value (for reference):
  // TODO: refactor this out

  double c = (double)(mlf->lpf[n-1]); // total number of factors
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
