#include <cinttypes>
#include <vector>
using namespace std;

#include "prelude.h"
#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

// input: fenwick tree array, size, index to increase, amount
void fw_add(vector<int64_t> &a, size_t idx, int64_t delta) {
  for (; idx < a.size(); idx |= idx + 1)
    a[idx] += delta;
}

// input: fenwick tree array, index, variable to sum into
int64_t fw_sum(vector<int64_t> &a, int64_t idx) {
  int64_t sum = 0;
  while (idx >= 0) {
    sum += a[idx];
    idx &= idx + 1;
    idx--;
  }
  return sum;
}

// input: fenwick tree array, inclusive left and right border
int64_t fw_from_to(vector<int64_t> &a, size_t l, size_t r) {
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
void mlComplexity(size_t w, size_t k, vector<double> &y, Fact &mlf, double gc) {
  size_t n = mlf.strLen;
  size_t entries = (n - w) / k + 1;

  //calculations (per nucleotide)
  double cMin = 2.0/n;           // at least 2 factors an any sequence, like AAAAAA.A
  double esl = expShulen(gc, n); // some wildly advanced estimation for avg. factor length
  double cAvg = 1.0 / (esl-1.0); // expected number of match length factors

  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = MIN(n, l + w) - 1;
    double cObs = (double)sumFromTo(mlf.lpf, l, r) / (double)w;
    y[j] = (cObs /* - cMin */) / (cAvg - cMin);
  }
}

// TODO: neue berechnung
void runComplexity(size_t w, size_t k, vector<double> &y, size_t n,
                   vector<list<Periodicity>> &ls) {
  size_t entries = (n - w) / k + 1;

  // for each position, get number of periodicities it is part of
  vector<int64_t> ft(n, 0);
  for (size_t i = 0; i < n; i++) {
    size_t len = 0;
    for (auto it = ls[i].begin(); it != ls[i].end(); it++) {
      Periodicity p = *it;
      size_t num = persFromRun(p); //# of period. corresponding to this run
      len += num;
      if (p.e + 1 < n)
        fw_add(ft, p.e + 1, -num);
    }
    fw_add(ft, i, len); // all in same list have same beginning
  }

  // collect result for each position from fenwick tree, sum up
  int64_t ppnMax = 0; // maximum observed periodicities per nucleotide
  int64_t *ps = new int64_t[n];
  for (size_t i = 0; i < n; i++) {
    ps[i] = fw_sum(ft, i);
    if (ps[i] > ppnMax)
      ppnMax = ps[i];
    if (i)
      ps[i] += ps[i - 1];
  }

  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = MIN(n, l + w) - 1;
    double ppnObs = (double)sumFromTo(ps, l, r) / w;
    y[j] = ppnObs / ppnMax;
  }
  delete[] ps;
}
