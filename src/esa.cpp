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
#include "prelude.h"

#include <cstdio>
#include <divsufsort64.h>

#include "esa.h"
#include "rmq.h"

// calculate suffix array using divsufsort
saidx_t *getSa(char const *seq, size_t n) {
  sauchar_t *t = (sauchar_t *)seq;
  saidx_t *sa = new saidx_t[n + 1];
  if (divsufsort(t, sa, (saidx_t)n) != 0) {
    printf("ERROR[esa]: suffix sorting failed.\n");
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
  saidx_t *sa = esa.sa;

  saidx_t *rank = new saidx_t[n]; // isa
  for (size_t i = 0; i < n; i++)
    rank[sa[i]] = i;
  esa.isa = rank;

  int64_t *lcp = new int64_t[n + 1];
  int64_t h = 0, j = 0;
  lcp[0] = lcp[n] = -1;
  for (size_t i = 0; i < n; i++) {
    if (rank[i] > 0) {
      j = sa[rank[i] - 1];
      while (t[i + h] == t[j + h]) {
        h++;
      }
      lcp[rank[i]] = h;
      if (h > 0)
        h--;
    }
  }
  esa.lcp = lcp;
}

Esa::Esa(char const *seq, size_t len) : str(seq), n(len) {
  this->sa = getSa(seq, n);
  calcLcp(*this);
}

Esa::~Esa() {
  delete[] this->sa;
  delete[] this->isa;
  delete[] this->lcp;
}

void Esa::print() const {
  printf("i\tSA\tLCP\tSuffix\n");
  for (size_t i = 0; i < this->n; i++)
    printf("%zu\t%d\t%ld\t%s\n", i, this->sa[i], this->lcp[i], this->str + this->sa[i]);
  printf("\t\t%ld\n", this->lcp[this->n]);
}

RMQ Esa::precomputeLcp() const { return RMQ(this->lcp, this->n + 1); }

int64_t Esa::getLcp(const RMQ &rmq, size_t sai, size_t saj) const {
  if (sai == saj)
    return this->n - sai;
  size_t l = MIN(this->isa[sai], this->isa[saj]) + 1;
  size_t r = MAX(this->isa[sai], this->isa[saj]);
  return rmq.get(l, r);
}
