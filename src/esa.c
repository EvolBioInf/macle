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

#include <stdio.h>
#include <divsufsort.h>

#include "eprintf.h"
#include "esa.h"
#include "rmq.h"

// calculate suffix array using divsufsort
int64_t *getSa(char *seq, size_t n) {
  sauchar_t *t = (sauchar_t *)seq;
  saidx_t *sa1 = (saidx_t *)emalloc(n * sizeof(saidx_t));
  if (divsufsort(t, sa1, (saidx_t)n) != 0) {
    printf("ERROR[esa]: suffix sorting failed.\n");
    exit(-1);
  }

  int64_t *sa2 = (int64_t *)emalloc(n * sizeof(int64_t));
  for (size_t i = 0; i < n; i++)
    sa2[i] = (long)sa1[i];
  free(sa1);
  return sa2;
}

/* calcLcp: compute LCP array using the algorithm in Figure 3
 *   of Kasai et al (2001). Linear-time longest-common-prefix
 *   computation in suffix arrays and its applications. LNCS 2089
 *   p. 191-192.
 */
void calcLcp(Esa *esa) {
  char *t = esa->str;
  size_t n = esa->n;
  int64_t *sa = esa->sa;

  int64_t *rank = (int64_t *)emalloc(n * sizeof(int64_t)); // isa
  for (size_t i = 0; i < n; i++)
    rank[sa[i]] = i;
  esa->isa = rank;

  int64_t *lcp = (int64_t *)emalloc((n + 1) * sizeof(int64_t));
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
  esa->lcp = lcp;
}

Esa *getEsa(char *seq, size_t n) {
  Esa *esa = (Esa *)emalloc(sizeof(Esa));
  esa->sa = getSa(seq, n);
  esa->str = seq;
  esa->n = n;
  calcLcp(esa);
  return esa;
}

void freeEsa(Esa *esa) {
  free(esa->sa);
  free(esa->isa);
  free(esa->lcp);
  free(esa);
}

void printEsa(Esa *esa) {
  printf("i\tSA\tLCP\tSuffix\n");
  for (size_t i = 0; i < esa->n; i++)
    printf("%zu\t%zu\t%ld\t%s\n", i, esa->sa[i], esa->lcp[i], esa->str + esa->sa[i]);
  printf("\t\t%ld\n", esa->lcp[esa->n]);
}

int64_t *precomputeLcp(Esa *esa) {
  return precomputeRMQ(esa->lcp, esa->n+1);
}

int64_t getLcp(Esa *esa, int64_t *lcptab, size_t sai, size_t saj) {
  if (sai==saj)
    return esa->n-sai;
  size_t l = MIN(esa->isa[sai], esa->isa[saj])+1;
  size_t r = MAX(esa->isa[sai], esa->isa[saj]);
  return RMQ(esa->lcp, esa->n+1, lcptab, l, r);
}
