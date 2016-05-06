/***** esa.c **************************************
 * Description: Enhanced Suffix Array.
 * Reference: Abouelhoda, Kurtz, and Ohlebusch
 *   (2002). The enhanced suffix array and its
 *   applications to genome analysis. Proceedings
 *   of the Second Workshop on Algorithms in
 *   Bioinformatics, Springer Verlag, Lectore Notes
 *   in Compter Science.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Mon Jul 15 11:11:19 2013
 **************************************************/
#include "bench.h"
#include <cinttypes>
#include <iostream>
#include <algorithm>
#include <vector>
using namespace std;

#include <divsufsort64.h>

#include "esa.h"
#include "rmq.h"

// calculate suffix array using divsufsort
vector<saidx64_t> getSa(char const *seq, size_t n) {
  sauchar_t *t = (sauchar_t *)seq;
  vector<saidx64_t> sa(n + 1);
  if (divsufsort64(t, sa.data(), (saidx64_t)n) != 0) {
    cout << "ERROR[esa]: suffix sorting failed." << endl;
    exit(-1);
  }
  return sa;
}

/* calcLcp: compute LCP array using the algorithm in Figure 3
 *   of Kasai et al (2001). Linear-time longest-common-prefix
 *   computation in suffix arrays and its applications. LNCS 2089
 *   p. 191-192.
 */
void calcLcp(Esa &esa) {
  char const *t = esa.str;
  size_t n = esa.n;
  auto &sa = esa.sa;
  auto &rank = esa.isa;

  esa.lcp = vector<saidx64_t>(n + 1);
  auto &lcp = esa.lcp;
  int64_t h = 0, j = 0;
  lcp[0] = lcp[n] = -1;
  for (size_t i = 0; i < n; i++) {
    if (rank[i] > 0) {
      j = sa[rank[i] - 1];
      while (t[i + h] == t[j + h])
        h++;
      lcp[rank[i]] = h;
      if (h > 0)
        h--;
    }
  }
}

// alternative lcp calculation (PLCP algorithm, Ohlebusch book)
// from: EvolBioInf/andi source code
void esa_init_LCP(Esa &esa) {
  const char *S = esa.str;
  auto &SA = esa.sa;
  saidx64_t len = esa.n;

  // The LCP array is one element longer than S.
  esa.lcp = vector<saidx64_t>(len + 1);
  auto &LCP = esa.lcp;

  LCP[0] = -1;
  LCP[len] = -1;

  // Allocate temporary arrays
  vector<saidx64_t> PHI(len);
  auto &PLCP = PHI;

  PHI[SA[0]] = -1;
  saidx64_t k;
  ssize_t i;

  for (i = 1; i < len; i++) {
    PHI[SA[i]] = SA[i - 1];
  }

  ssize_t l = 0;
  for (i = 0; i < len; i++) {
    k = PHI[i];
    if (k != -1) {
      while (S[k + l] == S[i + l]) {
        l++;
      }
      PLCP[i] = l;
      l--;
      if (l < 0)
        l = 0;
    } else {
      PLCP[i] = -1;
    }
  }

  // unpermutate the LCP array
  for (i = 1; i < len; i++)
    LCP[i] = PLCP[SA[i]];
}

Esa::Esa(char const *seq, size_t len) : str(seq), n(len) {
  tick();
  sa = getSa(seq, n);
  tock("libdivsufsort");

  isa = vector<saidx64_t>(n);
  for (size_t i = 0; i < n; i++)
    isa[sa[i]] = i;

  tick();
  // calcLcp(*this);
  esa_init_LCP(*this);
  tock("calcLCP");
}

void Esa::print() const {
  cout << "i\tSA\tISA\tLCP\tSuffix" << endl;
  for (size_t i = 0; i < this->n; i++)
    cout << i << "\t" << sa[i] << "\t" << isa[i] << "\t" << lcp[i] << "\t" << str+sa[i] << endl;
}

RMQ Esa::precomputeLcp() const { return RMQ(this->lcp); }

int64_t Esa::getLcp(const RMQ &rmq, size_t sai, size_t saj) const {
  if (sai == saj)
    return this->n - sai;
  size_t l = min(this->isa[sai], this->isa[saj]) + 1;
  size_t r = max(this->isa[sai], this->isa[saj]);
  return rmq.get(l, r);
}

//reduce esa to half (seq+$+revseq+$ -> seq+$) without recomputing
//important! asserting that n and str are replaced by user!
void reduceEsa(Esa &esa) {
  size_t n = esa.sa.size()/2;
  vector<saidx64_t> sa(n);
  vector<saidx64_t> isa(n);
  vector<int64_t> lcp(n+1);
  lcp[0]=lcp[n]=-1;
  sa[0]=n-1;
  isa[n-1]=0;
  int64_t lcptmp=INT64_MAX;
  size_t ind=1;
  for (size_t i=1; i<n*2; i++) {
    if ((size_t)esa.sa[i]<n-1) {
      sa[ind] = esa.sa[i];
      isa[sa[ind]] = ind;
      lcp[ind] = min(lcptmp, esa.lcp[i]);
      ind++;
      lcptmp=INT64_MAX;
    } else {
      lcptmp = min(lcptmp, esa.lcp[i]);
    }
  }
  esa.sa.clear();
  esa.isa.clear();
  esa.lcp.clear();
  esa.sa = sa;
  esa.isa = isa;
  esa.lcp = lcp;
}

