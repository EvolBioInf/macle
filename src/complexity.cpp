#include "eprintf.h"
#include "prelude.h"
#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

// input: fenwick tree array, size, index to increase, amount
void fw_add(int64_t *a, size_t sz, size_t idx, int64_t delta) {
  for (; idx < sz; idx |= idx + 1)
    a[idx] += delta;
}

// input: fenwick tree array, index, variable to sum into
int64_t fw_sum(int64_t *a, int64_t idx) {
  int64_t sum = 0;
  while (idx >= 0) {
    sum += a[idx];
    idx &= idx + 1;
    idx--;
  }
  return sum;
}

// input: fenwick tree array, inclusive left and right border
int64_t fw_from_to(int64_t *a, size_t l, size_t r) {
  return fw_sum(a, r) - (l ? fw_sum(a, l - 1) : 0);
}

// input: alphabet size, string length
// output: upper bound (not exact) on number of ML factors
uint64_t maxFacts(uint64_t g, uint64_t n) {
  uint64_t fact = 0; // counted possible ML Factors
  uint64_t sum = 0;  // current length of string
  uint64_t l = 1;    // current factor length
  while (1) {
    uint64_t num = pow(g, l);
    uint64_t result = sum + l * num;
    if (result > n)
      break;
    sum = result;
    fact += num;
    l++;
  }
  if (sum < n)
    fact += (n - sum) / l; // estimate remaining length left to fill
  return fact;
}

// calculate match length complexity for sliding windows
// input: sane w and k, allocated array for results, match length factors, gc content
void mlComplexity(size_t w, size_t k, double *y, Fact *mlf, double gc) {
  size_t n = mlf->strLen;
  size_t entries = (n - w) / k + 1;

  // TODO: does this calculation make sense? (probably not)
  double cMin = 2.0 / (double)w; // at least 2 factors an any sequence, like AAAAAA.A
  double esl = expShulen(gc, w); // some wildly advanced estimation for avg. factor length
  double cMax = 1. / (esl - 1.);
  /* double cMax = maxFacts(4, w); // worst case upper bound */

  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = MIN(n, l + w) - 1;
    double cObs = sumFromTo(mlf->lpf, l, r);
    y[j] = (cObs/w /*-cMin*/) / (cMax - cMin);
    /* y[j] = cObs / cMax; */
  }
}

// TODO: better ideas for analysis of runs?
/*
void runComplexity(size_t w, size_t k, double *y, size_t n, List **ls) {
  size_t entries = (n - w) / k + 1;

  // for each start point, calculate longest periodicity
  size_t *maxpl = ecalloc(n, sizeof(size_t));
  size_t absmax = 0;
  for (size_t i = 0; i < n; i++) {
    List *curr = ls[i];
    size_t currmax = 0;
    if (curr)
      do {
        Periodicity *p = (Periodicity *)curr->value;
        size_t pl = p->e - p->b + 1;
        if (pl > currmax)
          currmax = pl;
        curr = curr->next;
      } while (curr);
    if (currmax > absmax)
      absmax = currmax;
    maxpl[i] = currmax;
  }

  // double rMax = 1.029*w;
  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = MIN(n, l + w) - 1;

    // double rObs = sumFromTo(lens, l, r);
    // y[j] = rObs / rMax;

    size_t rMax = 0;
    for (size_t i = l; i <= r; i++)
      if (maxpl[i] > rMax)
        rMax = maxpl[i];
    y[j] = (double)rMax / absmax;
  }

  free(maxpl);
}
*/

void runComplexity(size_t w, size_t k, double *y, size_t n, List **ls) {
  size_t entries = (n - w) / k + 1;

  // for each position, get number of periodicities it is part of
  int64_t *ft = (int64_t*)ecalloc(n, sizeof(int64_t));
  for (size_t i = 0; i < n; i++) {
    size_t len = 0;
    for (eachListItem(curr, ls[i])) {
      Periodicity *p = (Periodicity *)curr->value;
      size_t num = persFromRun(p); //# of period. corresponding to this run
      len += num;
      if (p->e + 1 < n)
        fw_add(ft, n, p->e + 1, -num);
    }
    fw_add(ft, n, i, len); // all in same list have same beginning
  }

  // collect result for each position from fenwick tree, sum up
  int64_t ppnMax = 0; // maximum observed periodicities per nucleotide
  int64_t *ps = (int64_t*)ecalloc(n, sizeof(int64_t));
  for (size_t i = 0; i < n; i++) {
    ps[i] = fw_sum(ft, i);
    if (ps[i] > ppnMax)
      ppnMax = ps[i];
    if (i)
      ps[i] += ps[i - 1];
  }
  free(ft);

  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = MIN(n, l + w) - 1;
    double ppnObs = (double)sumFromTo(ps, l, r) / w;
    y[j] = ppnObs / ppnMax;
  }

  free(ps);
}
