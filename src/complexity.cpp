#include <queue>
#include <cassert>
#include <iostream>
#include <algorithm>
using namespace std;

#include "args.h"  //args.p
#include "bench.h" //ticktock
#include "complexity.h"
#include "index.h"
#include "matchlength.h"
#include "shulen.h"

// input: prefix-sum array, left and right bound (inclusive)
#define sumFromTo(a, l, r) ((a)[(r)] - ((l) ? (a)[(l)-1] : 0))

// number of sliding windows with width w that fit into n with step k
size_t numEntries(size_t n, size_t w, size_t k) {
  assert(w <= n);
  assert(k <= w);
  if (w==0) //1 entry in global mode (no window set)
    return 1;
  return (n - w) / k + 1;
}

// usage: for each_window(total, windowsize, step) -> sets j,l,r in each iteration
#define each_window(n, w, k)                                                             \
  (size_t numj = numEntries(n, w, k), j = 0, l = j * k, r = min(n, l + w) - 1; j < numj; \
   j++, l = j * k, r = min(n, l + w) - 1)

// given sequence length, desired window and step size and a list of intervals
// with invalid nucleotides, returns a list of windows (iteration numbers) that
// should be ignored in the complexity calculation
queue<size_t> calcNAWindows(size_t offset, size_t n, size_t w, size_t k,
                             vector<pair<size_t, size_t>> const &badiv) {
  queue<size_t> na;
  list<pair<size_t,size_t>> bad;
  auto currbad = badiv.begin();

  for
    each_window(n, w, k) {
      //kick out intervals that are now outside of window
      for (auto it=bad.begin(); it != bad.end(); it++)
        if (it->second < l+offset) {
          it = bad.erase(it);
          it--;
        }
      //add possible new bad intervals
      while (currbad != badiv.end() && currbad->first <= r+offset) {
        bad.push_back(*currbad);
        currbad++;
      }

      size_t sum = 0;
      for (auto &i : bad)
        sum += max(i.first, l+offset) - min(r+offset, i.second) + 1;
      if ((double)sum / (double)w > 0.05)
        na.push(j);
    }
  return na;
}

//get number of bad nucleotides in given interval of given sequence data
pair<size_t,size_t> numBad(size_t offset, size_t len, ComplexityData const &dat) {
  if (offset==0 && len==dat.len)
    return make_pair(dat.numbad, dat.bad.size()); //stored in data

  size_t sum = 0;
  size_t ivs = 0;
  auto it = lower_bound(dat.bad.begin(),dat.bad.end(),make_pair(offset,offset),
      [](pair<size_t,size_t> a, pair<size_t,size_t> b){return a.second < b.second;});
  if (it != dat.bad.begin() && it->first < offset)
    it--;

  while (it != dat.bad.end() && it->first < offset+len) {
    int64_t add = max(0L, (int64_t)min(offset+len, it->second) - (int64_t)max(it->first, offset) + 1L);
    // cerr << it->first << " - " << it->second << " -> " << add << endl;
    sum += add;
    ivs++;
    it++;
  }
  return make_pair(sum,ivs);
}

// calculate match length complexity for sliding windows
// input: sequence length, sane w and k, allocated array for results, extracted data
void mlComplexity(size_t offset, size_t n, size_t w, size_t k, vector<double> &y, ComplexityData const &dat) {
  bool globalMode = n==w;
  auto badpart = numBad(offset, n, dat);
  size_t numbad = badpart.first;
  size_t badivs = badpart.second;

  // compute observed number of match factors for every prefix
  vector<size_t> ps(n);
  size_t nextfact = lower_bound(dat.mlf.begin(),dat.mlf.end(),offset+1,[](size_t a,size_t b){return a<b;})-dat.mlf.begin();
  ps[0] = 1;
  for (size_t i = 1; i < n; i++) {
    ps[i] = ps[i - 1];
    if (nextfact < dat.mlf.size() && i+offset == dat.mlf[nextfact]) {
      ps[i]++;
      nextfact++;
    }
  }

  // calculations (per nucleotide)
  double cMin = 2.0 / (dat.len - dat.numbad); // at least 2 factors an any sequence, like AAAAAA.A

  // some wildly advanced estimation for avg. shulen length,
  // 2n because matches are from both strands
  double esl = expShulen(dat.gc, 2 * (dat.len - dat.numbad));
  if (globalMode) { //global complexity -> ignore NNNN... blocks, as if they are not there
    double fracbad = (double)numbad / (double)n;
    // cerr << fracbad << endl;
    if (fracbad > 0.05) {
      cerr << "WARNING: only " << (1-fracbad)*100 << "\% of sequence are valid DNA! "
           << "Ignoring " << badivs << " bad intervals..." << endl;
    }
  }

  // expected # of match length factors / nucleotide
  double cAvg = 1.0 / (esl - 1.0);
  // only subtract cMin if its not a degenerate case
  double cNorm = cAvg - cMin > 0 ? cAvg - cMin : cAvg;

  if (args.p) {
    cout  << "expected match factor length: " << esl-1 << endl;
    cout << "expected match factors per nucleotide: " << cNorm << endl;
  }

  // get bad window indices for given parameters
  queue<size_t> badj = calcNAWindows(offset, n, w, k, dat.bad);
  for
    each_window(n, w, k) {
      if (!globalMode)
        if (!badj.empty() && j == badj.front()) {
          y[j] = -1;
          badj.pop();
          continue;
        }

      size_t numfacs = sumFromTo(ps, l, r);
      double effectiveW = w;
      if (globalMode) {
        numfacs -= badivs;     //subtract number of bad intervals from total
        effectiveW -= numbad;  //we ignore the N-blocks
      }

      double cObs = (double)numfacs / effectiveW;
      y[j] = (cObs /* - cMin */) / cNorm;

      if (args.p) {
        cout << "observed match factors: " << numfacs << endl;
        cout << "observed match factors per nucleotide: " << cObs << endl;
        cout << "mlComplexity = avgPerNucl / estimated = "
              << cObs << "/" << cNorm << " = " << y[j] << endl;
      }
    }
}

ResultMat calcComplexities(size_t &w, size_t &k, Task task, ComplexityData const &dat) {
  bool globalMode = w==0;   // output one number (window = whole sequence)?
  bool wholeSeq = task.idx < 0;

  size_t offset = 0;
  size_t len = dat.len;
  if (!wholeSeq) {
    offset = dat.regions[task.idx].first;
    len    = dat.regions[task.idx].second;
  }
  if (task.end != task.start) {
    offset += task.start;
    len = task.end - task.start + 1;
  }

  // adapt window size and interval
  if (globalMode)
    w = k = len;
  w = min(w, len); // biggest window = whole seq.
  if (w != 0 && k == 0)
    k = max((size_t)1, w / 10); // default interval = 1/10 of window
  k = min(k, w);                // biggest interval = window size

  // cerr << offset << " " << len << " " << w << " " << k << endl;

  // array for results for all sequences in file
  size_t entries = numEntries(len, w, k);
  ResultMat ys(1, make_pair("", vector<double>(entries)));
  string name = wholeSeq ? dat.name :  dat.labels[task.idx];
  ys[0].first = name + " (MC)";

  tick();
  mlComplexity(offset, len, w, k, ys[0].second, dat);
  tock("mlComplexity");
  return ys;
}
