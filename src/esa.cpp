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
#include <cinttypes>
#include <iostream>
#include <algorithm>
using namespace std;

#include <divsufsort64.h>

#include "bench.h"
#include "esa.h"

using namespace sdsl;

// calculate suffix array using divsufsort
int_vector<VECBIT> getSa(char const *seq, size_t n) {
  sauchar_t *t = (sauchar_t *)seq;
  vector<saidx64_t> sa(n + 1);
  if (divsufsort64(t, sa.data(), (saidx64_t)n) != 0) {
    cout << "ERROR[esa]: suffix sorting failed." << endl;
    exit(-1);
  }
  int_vector<VECBIT> ret(n+1);
  for (size_t i=0; i<n+1; i++)
    ret[i] = sa[i];
  if (!VECBIT)
    util::bit_compress(ret);
  return ret;
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

  esa.lcp.resize(n + 1);
  auto &lcp = esa.lcp;
  int64_t h = 0, j = 0;
  lcp[0] = lcp[n] = 0;
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
  if (!VECBIT)
    util::bit_compress(lcp);
}

Esa::Esa(char const *seq, size_t len) : str(seq), n(len) {
  // string s(str);
  // construct_im(sa, s.c_str(), 1);
  tick();
  sa = getSa(seq, n);
  tock("libdivsufsort");

  isa = int_vector<VECBIT>(n+1);
  for (size_t i = 0; i < n; i++)
    isa[sa[i]] = i;
  if (!VECBIT)
    util::bit_compress(isa);

  tick();
  calcLcp(*this);
  tock("calcLCP");
}

void Esa::print() const {
  cout << "i\tSA\tISA\tLCP\tSuffix" << endl;
  for (size_t i = 0; i < this->n; i++)
    cout << i << "\t" << sa[i] << "\t" << isa[i] << "\t" << lcp[i] << "\t" << str + sa[i]
         << endl;
}

RMQ Esa::precomputeLcp() const { return RMQ(lcp); }

int64_t Esa::getLcp(RMQ const &rmq, size_t sai, size_t saj) const {
  if (sai == saj)
    return this->n - sai;
  size_t l = min(this->isa[sai], this->isa[saj]) + 1;
  size_t r = max(this->isa[sai], this->isa[saj]);
  return rmq(l, r);
}

// reduce esa to half (seq+$+revseq+$ -> seq+$) without recomputing
// important! asserting that n and str are replaced by user!
void reduceEsa(Esa &esa) {
  size_t n = esa.sa.size() / 2;
  int_vector<VECBIT> sa(n);
  int_vector<VECBIT> isa(n);
  int_vector<VECBIT> lcp(n + 1);
  lcp[0] = lcp[n] = 0;
  sa[0] = n - 1;
  isa[n - 1] = 0;
  int64_t lcptmp = INT64_MAX;
  size_t ind = 1;
  for (size_t i = 1; i < n * 2; i++) {
    if ((size_t)esa.sa[i] < n - 1) {
      sa[ind] = esa.sa[i];
      isa[sa[ind]] = ind;
      lcp[ind] = min(lcptmp, (int64_t)esa.lcp[i]);
      ind++;
      lcptmp = INT64_MAX;
    } else {
      lcptmp = min(lcptmp, (int64_t)esa.lcp[i]);
    }
  }
  esa.sa.resize(0);
  esa.isa.resize(0);
  esa.lcp.resize(0);
  if (!VECBIT) {
    util::bit_compress(sa);
    util::bit_compress(isa);
    util::bit_compress(lcp);
  }
  esa.sa = sa;
  esa.isa = isa;
  esa.lcp = lcp;
}
