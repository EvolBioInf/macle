/***** factor.c ***********************************
 * Description: Compute the longest previous factor
 *   array using a suffix array and a longest
 *   common prefix array.
 * Reference: Crochemore, M., Ilie, L. and Smyth,
 *   W. F. (2008). A simple algorithm for com-
 *   puting the Lempel Ziv factorization. In:
 *   Data Compression Conference, p. 482-488.
 *   Computing longest previous factor in linear
 *   time and applications.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 10:29:09 2013
 **************************************************/
#include "lempelziv.h"

#include <algorithm>
#include <stack>
#include <vector>
using namespace std;

// sum type used in lpf algorithm
struct Elem {
  Elem(int64_t i, int64_t lcp);
  int64_t i;
  int64_t lcp;
};
Elem::Elem(int64_t a, int64_t b) : i(a), lcp(b) {}

/*
 * computeLpf: Compute longest previous factor
 * Reference: M. Crochemore, L. Ilie, W.F. Smyth.
 *   A simple algorithm for computing the Lempel-Ziv
 *   factorization, in: J.A. Storer, M.W. Marcellini
 *   (Eds.), 18th Data Compression Conference, IEEE
 *   Computer Society, Los Alamitos, CA, 2008,
 *   pp. 482-488.
 */
void computeLpf(vector<size_t> &lpf, vector<int64_t> &prevOcc, Esa const &esa) {
  size_t n = esa.n;
  auto &sa = esa.sa;
  auto &lcp = esa.lcp;

  lpf = vector<size_t>(n + 1);
  prevOcc = vector<int64_t>(n);

  lpf[n] = 0;

  stack<Elem> s;
  s.push(Elem(0, 0));

  for (size_t i = 1; i <= n; i++) {
    int64_t currLcp = i == n ? 0 : lcp[i];
    int64_t sai = i == n ? -1 : sa[i];
    while (!s.empty() && sai < sa[s.top().i]) {
      Elem v = s.top();
      s.pop();
      lpf[sa[v.i]] = max(v.lcp, currLcp);
      currLcp = min(v.lcp, currLcp);
      // fill prevOcc
      if (lpf[sa[v.i]] == 0)
        prevOcc[sa[v.i]] = -1;
      else if (v.lcp > currLcp)
        prevOcc[sa[v.i]] = sa[s.top().i];
      else
        prevOcc[sa[v.i]] = sai;
    }
    if (i < n)
      s.push(Elem((int64_t)i, currLcp));
  }
}

// alternative prevOcc calculation (from same paper)
void computeLpf2(vector<size_t> &lpf, vector<int64_t> &prevOcc, Esa const &esa) {
  size_t n = esa.n;
  auto &sa = esa.sa;
  vector<int64_t> lprev(n);
  vector<int64_t> lnext(n);
  vector<int64_t> prevl(n);
  vector<int64_t> prevr(n);
  lpf = vector<size_t>(n, 0);
  vector<int64_t> lpfl(n, 0);
  vector<int64_t> lpfr(n, 0);
  prevOcc = vector<int64_t>(n);
  for (size_t i = 0; i < esa.n; i++) {
    lprev[sa[i]] = i == 0 ? -1 : sa[i - 1];
    lnext[sa[i]] = i == esa.n - 1 ? -1 : sa[i + 1];
  }
  for (int64_t sai = esa.n - 1; sai >= 0; sai--) {
    prevl[sai] = lprev[sai];
    prevr[sai] = lnext[sai];
    if (lprev[sai] != -1)
      lnext[lprev[sai]] = lnext[sai];
    if (lnext[sai] != -1)
      lprev[lnext[sai]] = lprev[sai];
  }

  prevOcc[0] = -1;
  for (size_t i = 1; i < esa.n; i++) {
    size_t j = max(lpfl[i - 1] - 1, (int64_t)0);
    size_t k = max(lpfr[i - 1] - 1, (int64_t)0);
    if (prevl[i] == -1)
      lpfl[i] = 0;
    else {
      while (esa.str[i + j] == esa.str[prevl[i] + j])
        j++;
      lpfl[i] = j;
    }
    if (prevr[i] == -1)
      lpfr[i] = 0;
    else {
      while (esa.str[i + k] == esa.str[prevr[i] + k])
        k++;
      lpfr[i] = k;
    }
    lpf[i] = max(lpfl[i], lpfr[i]);
    if (lpf[i] == 0)
      prevOcc[i] = -1;
    else if (lpfl[i] > lpfr[i])
      prevOcc[i] = prevl[i];
    else
      prevOcc[i] = prevr[i];
  }
}

void computeLZFact(Fact &lzf, Esa const &esa, bool alternative) {
  lzf.prevOcc.clear();
  lzf.fact.clear();
  lzf.lpf.clear();
  lzf.str = esa.str;
  lzf.strLen = esa.n;
  if (alternative)
    computeLpf2(lzf.lpf, lzf.prevOcc, esa);
  else
    computeLpf(lzf.lpf, lzf.prevOcc, esa);

  lzf.fact.push_back(0);
  while (lzf.fact.back() < esa.n)
    lzf.fact.push_back(lzf.fact.back() + max((size_t)1, lzf.lpf[lzf.fact.back()]));
  lzf.fact.pop_back();
}
