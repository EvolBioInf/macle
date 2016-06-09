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
void computeLpf(uint_vec &lpf, vector<int64_t> &prevOcc, Esa const &esa) {
  size_t n = esa.n;
  auto &sa = esa.sa;
  auto &lcp = esa.lcp;

  lpf = uint_vec(n + 1);
  prevOcc = vector<int64_t>(n);

  lpf[n] = 0;

  stack<Elem> s;
  s.push(Elem(0, 0));

  for (size_t i = 1; i <= n; i++) {
    int64_t currLcp = i == n ? 0 : (int64_t)lcp[i];
    int64_t sai = i == n ? -1 : (int64_t)sa[i];
    while (!s.empty() && sai < (int64_t)sa[s.top().i]) {
      Elem v = s.top();
      s.pop();
      lpf[(size_t)sa[v.i]] = max(v.lcp, currLcp);
      currLcp = min(v.lcp, currLcp);
      // fill prevOcc
      if (lpf[(size_t)sa[v.i]] == 0)
        prevOcc[(size_t)sa[v.i]] = -1;
      else if (v.lcp > currLcp)
        prevOcc[(size_t)sa[v.i]] = (int64_t)sa[s.top().i];
      else
        prevOcc[(size_t)sa[v.i]] = sai;
    }
    if (i < n)
      s.push(Elem((int64_t)i, currLcp));
  }
}

// alternative prevOcc calculation (from same paper)
void computeLpf2(uint_vec &lpf, vector<int64_t> &prevOcc, Esa const &esa) {
  size_t n = esa.n;
  auto &sa = esa.sa;
  vector<int64_t> lprev(n);
  vector<int64_t> lnext(n);
  vector<int64_t> prevl(n);
  vector<int64_t> prevr(n);
  lpf = uint_vec(n, 0);
  vector<int64_t> lpfl(n, 0);
  vector<int64_t> lpfr(n, 0);
  prevOcc = vector<int64_t>(n);
  for (size_t i = 0; i < n; i++) {
    lprev[(size_t)sa[i]] = i == 0 ? -1 : (int64_t)sa[i - 1];
    lnext[(size_t)sa[i]] = i == n - 1 ? -1 : (int64_t)sa[i + 1];
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
  lzf.lpf.resize(0);
  lzf.str = esa.str;
  lzf.strLen = esa.n;
  if (alternative)
    computeLpf2(lzf.lpf, lzf.prevOcc, esa);
  else
    computeLpf(lzf.lpf, lzf.prevOcc, esa);
#ifdef USE_SDSL
  sdsl::util::bit_compress(lzf.lpf);
#endif

  vector<size_t> lzftmp;
  lzftmp.push_back(0);
  while (lzftmp.back() < esa.n)
    lzftmp.push_back(lzftmp.back() + max(1UL, (size_t)lzf.lpf[lzftmp.back()]));
  lzftmp.pop_back();

  lzf.fact = uint_vec(lzftmp.size());
  for (size_t i=0; i<lzftmp.size(); i++)
    lzf.fact[i] = lzftmp[i];
#ifdef USE_SDSL
  sdsl::util::bit_compress(lzf.fact);
#endif
}
