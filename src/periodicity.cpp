#include "bench.h"

#include "lempelziv.h"
#include "periodicity.h"
#include "esa.h"

#include <algorithm>
#include <string>
#include <vector>
#include <list>
using namespace std;

Periodicity::Periodicity(size_t x, size_t y, size_t z) : b(x), e(y), l(z) {}
void printPeriodicity(Periodicity &p) { printf("(%zu,%zu,%zu)\n", p.b, p.e, p.l); }

// after some LCP is above this constant, the RMQ method is used
#define LAZY_CONST 10

// for small lcps naive is faster -> try naive, then fall back to efficient RMQ
int64_t lcp(const Esa &esa, const RMQ &tab, size_t i, size_t j) {
  i--;
  j--;
  if (i >= esa.n || j >= esa.n)
    return 0;
  uint64_t k = 0;
  bool lazy = false;
  while (!lazy && max(i, j) + k < esa.n && esa.str[k + i] == esa.str[k + j])
    if (++k > LAZY_CONST)
      lazy = true;
  if (!lazy)
    return k;
  k = esa.getLcp(tab, i, j);
  return max(k, (uint64_t)0);
}

int64_t lcs(const Esa &resa, const RMQ &rlcptab, size_t i, size_t j) {
  int64_t maxval = lcp(resa, rlcptab, resa.n - i + 1, resa.n - j + 1);
  return max(maxval, (int64_t)0);
}

// naive: length of lcp of suffixes i and j of given seq (1-indexed)
size_t lcp2(char const *str, size_t n, size_t i, size_t j) {
  i--;
  j--;
  uint64_t k = 0;
  while (max(i, j) + k < n && str[k + i] == str[k + j])
    k++;
  return k;
}

// naive: length of lcs of prefixes 1..i and 1..j of given seq (input 1-indexed)
size_t lcs2(char const *str, size_t i, size_t j) {
  i--;
  j--;
  int64_t k = 0;
  while (((int64_t)min(i, j)) - k >= 0 && str[i - k] == str[j - k])
    k++;
  return k;
}

// max. periodicities in O(nlog(n)) (Algorithm 5.16), mainly for reference / comparison
vector<Periodicity> getPeriodicities2(Esa const &esa) {
  size_t n = esa.n;
  char const *str = esa.str;
  vector<Periodicity> ps;
  ps.reserve(n);
  for (size_t l = 1; l <= n / 2; l++) {
    size_t i = 2 * l + 1;
    while (i <= n + 1) {
      size_t L = lcs2(str, i - 1 - l, i - 1);
      size_t R = lcp2(str, n, i - l, i);
      if (L + R >= l && R < l) {
        size_t b = i - l - L - 1; // again 0-indexed
        size_t e = i - 1 + R - 1;
        ps.push_back(Periodicity(b, e, l));
      }
      i += l;
    }
  }
  ps.shrink_to_fit();
  return ps;
}

// start position (1-indexed)
static inline size_t factStart(Fact const &f, size_t i) { return f.fact[i] + 1; }

// end position (1-indexed)
static inline size_t factEnd(Fact const &f, size_t i) {
  return (i < f.fact.size() - 1 ? f.fact[i + 1] : f.strLen);
}

// helper to avoid duplication
static inline bool addPeriodicity(bool runsOnly, PerLists &Lt1,
                                  Periodicity p) {
  if (!Lt1[p.b].empty() && runsOnly) { // if we want only runs,
    Periodicity last = Lt1[p.b].back();
    if (last.e != p.e) { // add only if prev. is not same interval
      Lt1[p.b].push_back(p);
      return true;
    }
  } else { // otherwise, add in any case
    Lt1[p.b].push_back(p);
    return true;
  }
  return false;
}

// Algorithm 5.17 - calculate type1 periodicities
PerLists calcType1Periodicities(bool runsOnly, Fact const &lzf,
                                                 Esa const &esa, size_t &pnum) {
  tick();
  auto lcptab = esa.precomputeLcp();
  tock("precomputeLcp");
  tick();
  string srev = esa.str;
  reverse(srev.begin(), srev.end());
  Esa revesa(srev.c_str(), srev.size());
  auto revlcptab = revesa.precomputeLcp();
  tock("prepare reverse seq");

  pnum = 0;
  // as required, array of lists of max. per. of type 1
  // indexed by start position, each sorted by end position
  PerLists Lt1(lzf.strLen + 1);
  for (size_t j = 1; j < lzf.fact.size(); j++) {
    size_t sjLen = factLen(lzf, j);       // |s_j|
    size_t sjm1Len = factLen(lzf, j - 1); // |s_(j-1)|
    size_t bj = factStart(lzf, j);        // b_j
    size_t bjm1 = factStart(lzf, j - 1);  // b_(j-1)
    size_t ej = factEnd(lzf, j);          // e_j
    size_t ejm1 = factEnd(lzf, j - 1);    // e_(j-1)

    size_t max = min(sjLen + sjm1Len - 1, ejm1);
    for (size_t l = 1; l <= max; l++) {
      /* size_t L = lcs2(lzf.str, ejm1 - l, ejm1); */
      size_t L = lcs(revesa, revlcptab, ejm1 - l, ejm1);
      /* size_t R = lcp2(lzf.str, lzf.strLen, bj - l, bj); */
      size_t R = lcp(esa, lcptab, bj - l, bj);
      if (L + R >= l && (R >= 1 || bj - l - L > bjm1)) {
        size_t b = bj - l - L - 1;
        size_t e = ejm1 + R - 1;
        if (addPeriodicity(runsOnly, Lt1, Periodicity(b, e, l)))
          pnum++;
      }
    }

    for (size_t l = 1; l <= sjLen; l++) {
      /* size_t L = lcs2(lzf->str, ejm1, ejm1 + l); */
      size_t L = lcs(revesa, revlcptab, ejm1, ejm1 + l);
      /* size_t R = lcp2(lzf->str, lzf->strLen, bj, bj + l); */
      size_t R = lcp(esa, lcptab, bj, bj + l);
      if (L + R >= l && bj + l - 1 + R <= ej && L < l) {
        size_t b = bj - L - 1;
        size_t e = ejm1 + l + R - 1;
        if (addPeriodicity(runsOnly, Lt1, Periodicity(b, e, l)))
          pnum++;
      }
    }
  }

  return Lt1;
}

// prevOcc conforming to algorithm: 1-indexed and current value instead of -1
static inline int64_t getPrevOcc(Fact const &lzf, size_t i) {
  int64_t tmp = lzf.prevOcc[lzf.fact[i]];
  return (tmp == -1 ? (int64_t)lzf.fact[i] : tmp) + 1;
}

// add given periodicity into corresponding same-beginning list
// so that list remains sorted by end positions. This hack is necessary!
// Ohlebusch/Algorithm 5.18 says "prepend". This is WRONG!!! Was a fun night
// debugging -.-'
void listInsert(PerLists &ls, Periodicity p) {
  auto &l = ls[p.b];
  if (l.empty() || l.front().e >= p.e) {
    l.insert(l.begin(), p);
    return;
  }
  auto it = l.begin();
  while (it != l.end() && (*it).e < p.e)
    ++it;
  l.insert(it, p);
}

// Algorithm 5.18 - calculate type 2 periodicities (proper substrings of LZ-factors)
void calcType2Periodicities(PerLists &Lt1, Fact const &lzf,
                            size_t &pnum) {
  for (size_t j = 1; j < lzf.fact.size(); j++) {
    size_t bj = factStart(lzf, j); // b_j
    size_t ej = factEnd(lzf, j);   // e_j
    if (factLen(lzf, j) >= 4) {
      int64_t dj = bj - getPrevOcc(lzf, j);
      for (size_t i = bj; i < ej - 1; i++) {
        for (auto p : Lt1[i - dj]) {
          if (p.e + dj >= ej - 1)
            break;
          /* listPrepend(&Lt1[i], newp); */ // XXX: doesn't work with prepend!!!
          // need to sort into list correctly!!!!!!!!!!!!!!!
          listInsert(Lt1, Periodicity(i, p.e + dj, p.l));
          pnum++;
        }
      }
    }
  }
}

// max. periodicities in O(n) using LZ-Factors, returns array of lists
PerLists getPeriodicityLists(bool runsOnly, Fact const &lzf,
                                              Esa const &esa, size_t &pnum) {
  auto Lt1 = calcType1Periodicities(runsOnly, lzf, esa, pnum);
  calcType2Periodicities(Lt1, lzf, pnum);
  return Lt1;
}

vector<Periodicity> collectPeriodicities(PerLists &pl, size_t pnum) {
  vector<Periodicity> ps;
  ps.reserve(pnum);
  for (size_t i = 0; i < pl.size(); i++) {
    for (auto p : pl[i])
      ps.push_back(p);
    pl[i].clear();
  }
  pl.clear();
  return ps;
}

vector<Periodicity> getPeriodicities(bool runsOnly, Fact const &lzf, Esa const &esa,
                                     size_t &pnum) {
  auto pl = getPeriodicityLists(runsOnly, lzf, esa, pnum);
  return collectPeriodicities(pl, pnum);
}
