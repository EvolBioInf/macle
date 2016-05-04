#include <cinttypes>
#include <vector>
#include <algorithm>
using namespace std;

#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

// calculate match length complexity for sliding windows
// input: sane w and k, allocated array for results, match length factors, gc content
void mlComplexity(size_t w, size_t k, vector<double> &y, Fact &mlf, double gc) {
  size_t n = mlf.strLen / 2; // mlf was calculated on both strands, we look on first only
  size_t entries = (n - w) / k + 1;

  // calculations (per nucleotide)
  double cMin = 2.0 / n; // at least 2 factors an any sequence, like AAAAAA.A
  // some wildly advanced estimation for avg. shulen length,
  // 2n because matches are from both strands
  double esl = expShulen(gc, 2 * n);

  double cAvg = 1.0 / (esl - 1.0); // expected # of match length factors / nucleotide

  // calc. for each window
  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = min(n, l + w) - 1;
    double cObs = (double)sumFromTo(mlf.lpf, l, r) / (double)w;
    y[j] = (cObs /* - cMin */) / (cAvg - cMin);
  }
}

void runComplexity(size_t w, size_t k, vector<double> &y, size_t n,
                   vector<list<Periodicity>> &ls) {
  size_t entries = (n - w) / k + 1;

  vector<int64_t> ps(n, 0);
  for (size_t i = 0; i < n; i++) {
    for (auto it = ls[i].begin(); it != ls[i].end(); it++) {
      Periodicity p = *it;
      if (perLen(p) < 4)
        continue;
      for (size_t j = p.b; j <= p.e; j += p.l) // increment for each starting "atom"
        ps[j]++;
    }
  }

  // prefix sum (for fast range sum retrieval)
  for (size_t i = 1; i < n; i++)
    ps[i] += ps[i - 1];

  double pMax = w; // no more atoms than window size
  for (size_t j = 0; j < entries; j++) {
    size_t l = j * k;
    size_t r = min(n, l + w) - 1;
    double pObs = sumFromTo(ps, l, r);
    y[j] = pObs / pMax;
  }
}
