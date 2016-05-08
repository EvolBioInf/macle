#include <cassert>
#include <cinttypes>
#include <vector>
#include <algorithm>
using namespace std;

#include "complexity.h"
#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"
#include "IntervalTree.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

size_t numEntries(size_t n, size_t w, size_t k) {
  assert(w <= n);
  assert(k > 0);
  return (n - w) / k + 1;
}

// usage: for each_window(total, windowsize, step) -> sets j,l,r in each iteration
#define each_window(n, w, k)                                                             \
  (size_t max = numEntries(n, w, k), j = 0, l = j * k, r = min(n, l + w) - 1; j < max;   \
   j++, l = j * k, r = min(n, l + w) - 1)

// calculate match length complexity for sliding windows
// input: sequence length, sane w and k, allocated array for results,
//        match length factors, gc content
void mlComplexity(size_t n, size_t w, size_t k, vector<double> &y,
                  vector<size_t> const &fact, double gc) {
  // compute observed number of match factors for every prefix
  vector<size_t> ps(n);
  size_t nextfact = 1;
  ps[0] = 1;
  for (size_t i = 1; i < n; i++) {
    ps[i] = ps[i - 1];
    if (nextfact < fact.size() && i == fact[nextfact]) {
      ps[i]++;
      nextfact++;
    }
  }

  // calculations (per nucleotide)
  double cMin = 2.0 / n; // at least 2 factors an any sequence, like AAAAAA.A
  // some wildly advanced estimation for avg. shulen length,
  // 2n because matches are from both strands
  double esl = expShulen(gc, 2 * n);
  double cAvg = 1.0 / (esl - 1.0); // expected # of match length factors / nucleotide
  // only subtract cMin if its not a degenerate case
  double cNorm = cAvg - cMin > 0 ? cAvg - cMin : cAvg;

  for
    each_window(n, w, k) {
      double cObs = (double)sumFromTo(ps, l, r) / (double)w;
      y[j] = (cObs /* - cMin */) / cNorm;
    }
}

void runComplexityOld(size_t n, size_t w, size_t k, vector<double> &y,
                      vector<list<Periodicity>> const &ls) {
  vector<int64_t> ps(n, 0);
  for (auto &l : ls) {
    for (auto p : l) {
      if (perLen(p) < FILTER)
        continue;
      for (size_t j = p.b; j <= p.e; j++) // mark nucleotides inside runs
        ps[j] |= 1;
    }
  }

  // prefix sum (for fast range sum retrieval)
  for (size_t i = 1; i < n; i++)
    ps[i] += ps[i - 1];

  for
    each_window(n, w, k) {
      // complexity = fraction of nucleotides outside
      // of any remaining runs (after filtering)
      double pObs = (double)sumFromTo(ps, l, r) / (double)w;
      y[j] = 1 - pObs;
    }
}

// get "information content" of window
void runComplexity(size_t n, size_t w, size_t k, vector<double> &y,
                   vector<list<Periodicity>> const &ls) {
  vector<int64_t> ps(n, 1); // all nucleotides marked
  vector<Interval<size_t>> intervals;
  for (auto &l : ls) {
    for (auto p : l) {
      if (perLen(p) < FILTER)
        continue;
      for (size_t j = p.b; j <= p.e; j++) // un-mark nucleotides inside runs
        ps[j] = 0;
      intervals.push_back(Interval<size_t>(p.b, p.e, p.l));
    }
  }

  // prefix sum over non-run nucleotides (for fast range sum retrieval)
  for (size_t i = 1; i < n; i++)
    ps[i] += ps[i - 1];
  // construct interval tree to quickly get runs overlapping the window
  IntervalTree<size_t> tree = IntervalTree<size_t>(intervals);

  // subtract 1 from w to cap result at 1, max to prevent div. by 0
  double pMax = max(1.0, (double)w - 1.0);
  for
    each_window(n, w, k) {
      size_t pObs = sumFromTo(ps, l, r); // count non-run nucl. in window
      // cout << j << ": " << pObs;
      vector<Interval<size_t>> runs;
      tree.findOverlapping(l, r, runs);
      for (auto &run : runs) {
        pObs += run.value; // add period lengths of runs touching window
        // cout << " + " << run.value;
      }
      // cout << " = " << pObs << endl;
      pObs = min(w, pObs); // we look at all overl. runs, could lead to sum > w
      // subtract 1 from pObs to make 0 possible (e.g. for AAAA..)
      y[j] = ((double)pObs - 1.0) / pMax;
    }
}
