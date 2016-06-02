#include <cassert>
#include <cmath>
#include <cinttypes>
#include <vector>
#include <algorithm>
#include <utility>
#include <queue>
using namespace std;

#include "args.h"  //args.p
#include "bench.h" //ticktock
#include "complexity.h"
#include "index.h"
#include "matchlength.h"
#include "periodicity.h"
#include "shulen.h"
#include "util.h"

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
      //kick out runs that are now outside of window
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
      [](pair<size_t,size_t> a, pair<size_t,size_t> b){return a.first < b.first;});
  if (it != dat.bad.begin() && it->first < offset)
    it--;

  while (it != dat.bad.end() && it->first < offset+len) {
    sum += min(offset+len, it->second) - max(it->first, offset) + 1;
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
  double gc = dat.gc;
  if (globalMode) { // adapt gc content (as if no N-blocks present)
    // oldgc = c(gc) / (c(n)+c(at)+c(gc))=seqlen -> newgc = oldgc * seqlen / (seqlen-c(n))
    gc = gc * dat.len / (dat.len-numbad);
  }

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
  double cMin = 2.0 / dat.len; // at least 2 factors an any sequence, like AAAAAA.A
  if (globalMode) // ignore N-blocks from sequence
    cMin = cMin * dat.len / (dat.len-numbad);

  // some wildly advanced estimation for avg. shulen length,
  // 2n because matches are from both strands
  double esl = 0;
  if (n != w)
    esl = expShulen(gc, 2 * dat.len);
  else { //global complexity -> ignore NNNN... blocks, as if they are not there
    esl = expShulen(gc, 2 * (dat.len - dat.numbad));

    double fracbad = (double)numbad / (double)n;
    if (globalMode && fracbad > 0.05) {
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
      //in global mode -> subtract number of bad intervals from total
      if (globalMode)
        numfacs -= badivs;

      //in global mode we ignore the N-blocks
      double effectiveW = n==w ? w-numbad : w;

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


// simulate complexity calculation on random dna sequence
// to get expected value. the value is invariant under seq. length
double calcAvgRunComplexity(size_t len, double gc, size_t reps) {
  double c = 0;
  vector<double> y(1);
  for (size_t rep=0; rep<reps; rep++) {
    string seq = randSeq(len, gc) + "$";
    ComplexityData dat;
    dat.pl = getRuns(seq);
    dat.gc = gc;
    runComplexity(0, len, len, len, y, dat, false);
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

void updateRunQueue(list<Periodicity> &runs, PerLists const &ls, size_t l, size_t r) {
  //kick out runs that are now outside of window
  for (auto it=runs.begin(); it != runs.end(); it++)
    if (it->e < l) {
      it = runs.erase(it);
      it--;
    }
  //add new runs
  size_t start = runs.empty() ? l : max(l, runs.back().b+1);
  for (size_t i=start; i<=r; i++)
    for (auto &p : ls[i])
      runs.push_back(p);
}

// get "information content" of window
void runComplexity(size_t offset, size_t n, size_t w, size_t k, vector<double> &y,
                   ComplexityData const &dat, bool calcAvg) {
  bool globalMode = n==w;
  auto badpart = numBad(offset, n, dat);
  size_t numbad = badpart.first;
  size_t badivs = badpart.second;
  double gc = dat.gc;
  if (globalMode) { // adapt gc content (as if no N-blocks present)
    // oldgc = c(gc) / (c(n)+c(at)+c(gc))=seqlen -> newgc = oldgc * seqlen / (seqlen-c(n))
    gc = gc * dat.len / (dat.len-numbad);
  }

  vector<int64_t> ps(n, 1); // all nucleotides marked
  for (size_t i=offset; i<offset+n; i++)
    for (auto &p : dat.pl[i])
      for (size_t j = p.b; j <= min(offset+n-1, p.e); j++) { // un-mark nucleotides inside runs
        ps[j-offset] = 0;
      }

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
  queue<size_t> badj = calcNAWindows(offset, n, w, k, dat.bad);
  for
    each_window(n, w, k) {
      if (!globalMode)
        if (!badj.empty() && j == badj.front()) {
          y[j] = -1;
          badj.pop();
          continue;
        }

      size_t info = sumFromTo(ps, l, r); // count non-run nucl. in window

      //in global mode, N-blocks are ignored, they have period length 1
      //-> subtract them away
      if (globalMode)
        info -= badivs;

      if (args.p)
        cout << "NuclNotInRuns: " << info << endl;

      // update list of runs that overlap the current window
      updateRunQueue(runs, dat.pl, l+offset, r+offset);

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

      //in global mode we ignore the N-blocks
      double effectiveW = n==w ? w - numbad : w;

      // we look at all overlapping runs, could lead to sum > w -> min(w,*)
      // subtract 1 from pObs to make 0 possible (e.g. for AAAA..)
      double pObs = (min(w, info) - 1.0)/effectiveW;
      y[j] = pObs / pAvg;

      if (args.p) {
        cout << "avgPerNucl = (infoContent - 1) / seqLen = " << pObs << endl;
        cout << "runComplexity = avgPerNucl / estimated = "
             << pObs << "/" << pAvg << " = " << y[j] << endl;
      }
    }
}

ResultMat calcComplexities(size_t &w, size_t &k, char m, size_t seqnum, vector<ComplexityData> const &dat) {
  bool globalMode = w==0;   // output one number (window = each whole sequence)?
  bool isJoined = dat[0].regions.size()>1; //which kind is the input data?
  size_t containedSeqs = max(dat[0].labels.size(), dat.size()); //how many (sub-)sequences are there?

  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  if (seqnum || isJoined)
    maxlen = minlen = seqnum && isJoined ? dat[0].regions[seqnum-1].second : dat[0].len;
  else //we need to consider the separate sequences, non-joined mode
    for (auto &d : dat) {
      maxlen = max(maxlen, d.len);
      minlen = min(minlen, d.len);
    }

  // adapt window size and interval
  w = min(w, minlen); // biggest window = whole (smallest) seq.
  if (w != 0 && k == 0)
    k = max((size_t)1, w / 10); // default interval = 1/10 of window
  k = min(k, w);                // biggest interval = window size

  // array for results for all sequences in file
  size_t entries = numEntries(maxlen, w, k);
  size_t numMetrics = m == 'b' ? 2 : 1; //complexity arrays per seq.
  size_t usedSeqs = seqnum || isJoined ? 1 : containedSeqs;
  ResultMat ys(numMetrics * usedSeqs, make_pair("", vector<double>(entries)));

  int col = 0;
  int start = seqnum ? max(seqnum-1,(size_t)0)       : 0;
  int end   = seqnum ? min(seqnum-1,containedSeqs-1) : (isJoined ? start : containedSeqs-1);
  for (int i=start; i<=end; i++) {
    auto &data    = isJoined ? dat[0]                                           : dat[i];
    string name   = isJoined ? (seqnum ? dat[0].labels[i] : dat[0].name)        : dat[i].name;
    size_t offset = isJoined ? (seqnum ? dat[0].regions[i].first : 0)           : 0;
    size_t len    = isJoined ? (seqnum ? dat[0].regions[i].second : dat[0].len) : dat[i].len;
    size_t currw  = globalMode ? len : w;
    size_t currk  = globalMode ? len : k;

    if (m != 'r') {
      tick();
      ys[col].first = name + " (MC)";
      mlComplexity(offset, len, currw, currk, ys[col].second, data);
      tock("mlComplexity");
      col++;
    }
    if (m != 'm') {
      tick();
      ys[col].first = name + " (RC)";
      runComplexity(offset, len, currw, currk, ys[col].second, data, true);
      tock("runComplexity");
      col++;
    }
  }
  return ys;
}
