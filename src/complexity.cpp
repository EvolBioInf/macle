#include <cassert>
#include <cmath>
#include <cinttypes>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>
using namespace std;

#include "args.h"
#include "complexity.h"
#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"
#include "util.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

size_t numEntries(size_t n, size_t w, size_t k) {
  assert(w <= n);
  assert(k > 0);
  return (n - w) / k + 1;
}

// usage: for each_window(total, windowsize, step) -> sets j,l,r in each iteration
#define each_window(n, w, k)                                                             \
  (size_t numj = numEntries(n, w, k), j = 0, l = j * k, r = min(n, l + w) - 1; j < numj; \
   j++, l = j * k, r = min(n, l + w) - 1)

// given sequence length, desired window and step size and a list of intervals
// with invalid nucleotides, returns a list of windows (iteration numbers) that
// should be ignored in the complexity calculation
queue<size_t> calcNAWindows(size_t n, size_t w, size_t k,
                             vector<pair<size_t, size_t>> const &badiv) {
  queue<size_t> na;
  list<pair<size_t,size_t>> bad;
  auto currbad = badiv.begin();

  for
    each_window(n, w, k) {
      //kick out runs that are now outside of window
      for (auto it=bad.begin(); it != bad.end(); it++)
        if (it->second < l) {
          it = bad.erase(it);
          it--;
        }
      //add possible new bad intervals
      while (currbad != badiv.end() && currbad->first <= r) {
        bad.push_back(*currbad);
        currbad++;
      }

      size_t sum = 0;
      for (auto &i : bad)
        sum += max(i.first, l) - min(r, i.second) + 1;
      if ((double)sum / (double)w > 0.05)
        na.push(j);
    }
  return na;
}

// calculate match length complexity for sliding windows
// input: sequence length, sane w and k, allocated array for results,
//        match length factors, gc content
void mlComplexity(size_t n, size_t w, size_t k, vector<double> &y,
                  vector<size_t> const &fact, double gc, vector<pair<size_t,size_t>> const &badiv) {
  //calculate number of bad nucleotides for global mode
  size_t numbad=0;
  if (n==w)
    for (auto &bad : badiv)
      numbad += bad.second - bad.first + 1;

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
  double esl = 0;
  if (n != w)
    esl = expShulen(gc, 2 * n);
  else { //global complexity -> ignore NNNN... blocks, as if they are not there
    esl = expShulen(gc, 2 * (n - numbad));

    double fracbad = (double)numbad / (double)n;
    if (n==w && fracbad > 0.05) {
      cerr << "WARNING: only " << (1-fracbad)*100 << "\% of sequence are valid DNA! "
           << "Ignoring " << badiv.size() << " bad intervals..." << endl;
    }
  }
  // expected # of match length factors / nucleotide
  double cAvg = 1.0 / (esl - 1.0);
  // only subtract cMin if its not a degenerate case
  double cNorm = cAvg - cMin > 0 ? cAvg - cMin : cAvg;

  // get bad window indices for given parameters
  queue<size_t> badj = calcNAWindows(n, w, k, badiv);
  for
    each_window(n, w, k) {
      if (n!=w)
        if (!badj.empty() && j == badj.front()) {
          y[j] = -1;
          badj.pop();
          continue;
        }

      size_t numfacs = sumFromTo(ps, l, r);
      //in global mode -> subtract number of bad intervals from total
      if (n==w)
        numfacs -= badiv.size();

      double cObs = (double)numfacs / (double)w;
      y[j] = (cObs /* - cMin */) / cNorm;
    }
}


// simulate complexity calculation on random dna sequence
// to get expected value. the value is invariant under seq. length
double calcAvgRunComplexity(size_t len, double gc, size_t reps) {
  double c = 0;
  vector<double> y(1);
  for (size_t rep=0; rep<reps; rep++) {
    string seq = randSeq(len, gc) + "$";
    auto ls = getRuns(seq);
    runComplexity(len, len, len, y, getRuns(seq), gc, vector<pair<size_t,size_t>>(0), false);
    c += y[0];
  }
  return c/reps;
}
// polynomial that was interpolated from data gathered via calcAvgRunComplexity
double estimateAvgRunComplexity(double gc) {
  // return 4.87*pow(gc, 4) - 9.768*pow(gc, 3) + 6.41*pow(gc, 2) - 1.516*gc + 0.766;
  return 12.363*pow(gc,6) - 37.167*pow(gc,5) + 47.21*pow(gc,4)
      - 32.397*pow(gc,3) + 12.064*pow(gc,2) - 2.075*gc + 0.778;
  // return 26.615*pow(gc,8) - 106.236*pow(gc,7) + 185.481*pow(gc,6) - 185.114*pow(gc,5)
  //     + 118.019*pow(gc,4) - 51.128*pow(gc,3) + 14.576*pow(gc,2) - 2.212*gc + 0.78;
}

void updateRunQueue(list<Periodicity> &runs, vector<list<Periodicity>> const &ls, size_t l, size_t r) {
  //kick out runs that are now outside of window
  for (auto it=runs.begin(); it != runs.end(); it++)
    if (it->e < l) {
      it = runs.erase(it);
      it--;
    }
  //add new runs
  for (size_t i=max(l,runs.back().b+1); i<=r; i++)
    for (auto &p : ls[i])
      runs.push_back(p);
}

// get "information content" of window
void runComplexity(size_t n, size_t w, size_t k, vector<double> &y,
                   vector<list<Periodicity>> const &ls, double gc,
                   vector<pair<size_t,size_t>> const &badiv, bool calcAvg) {
  vector<int64_t> ps(n, 1); // all nucleotides marked
  for (auto &l : ls)
    for (auto &p : l)
      for (size_t j = p.b; j <= p.e; j++) // un-mark nucleotides inside runs
        ps[j] = 0;

  // prefix sum over non-run nucleotides (for fast range sum retrieval)
  for (size_t i = 1; i < n; i++)
    ps[i] += ps[i - 1];

  // subtract 1 from w to cap result at 1, max to prevent div. by 0
  // divide by window length -> unnormalized max. information per nucleotide
  // this is only used directly for the expected value simulations
  double pAvg = max(1.0, (double)w - 1.0) / (double)w;
  //obtain expected nucleotid information (should be used as default!)
  if (calcAvg)
    pAvg = estimateAvgRunComplexity(gc);

  if (args.p) { //some user information
    cout << "GC-content: "<<gc<<endl;
    cout << "Estimated RC/nucleotide (via interpolated polynomial): "
      << "12.363*gc^6-37.167*gc^547.21*gc^4-32.397*gc^3+12.064*gc^2-2.075*gc+0.778 = "
      << pAvg << endl;
  }

  list<Periodicity> runs; //current runs overlapping window

  // get bad window indices for given parameters
  queue<size_t> badj = calcNAWindows(n, w, k, badiv);
  for
    each_window(n, w, k) {
      if (n!=w)
        if (!badj.empty() && j == badj.front()) {
          y[j] = -1;
          badj.pop();
          continue;
        }

      size_t info = sumFromTo(ps, l, r); // count non-run nucl. in window

      if (args.p)
        cout << "NuclNotInRuns: " << info << endl;

      // update list of runs that overlap the current window
      updateRunQueue(runs, ls, l, r);

      // add period lengths of runs touching window
      if (args.p)
        cout << "Add unique info (=period length) of runs: NuclNotInRuns";
      for (auto &run : runs) {
        info += run.l;
        if (args.p)
          cout << " + " << run.l;
      }
      if (args.p)
        cout << " = " << info << " = infoContent" << endl;

      // we look at all overlapping runs, could lead to sum > w -> min(w,*)
      // subtract 1 from pObs to make 0 possible (e.g. for AAAA..)
      double pObs = (min(w, info) - 1.0)/(double)w;

      if (args.p)
        cout << "avgPerNucl = (infoContent - 1) / seqLen = " << pObs << endl;

      y[j] = pObs / pAvg;

      if (args.p)
        cout << "runComplexity = avgPerNucl / estimated = "
             << pObs << "/" << pAvg << " = " << y[j] << endl;
    }
}
