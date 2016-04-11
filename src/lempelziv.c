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
#include "prelude.h"

#include "lempelziv.h"
#include "kvec.h"
#include "eprintf.h"

// sum type used in lpf algorithm
typedef struct elem {
  int64_t i;
  int64_t lcp;
} Elem;

/*
 * computeLpf: Compute longest previous factor
 * Reference: M. Crochemore, L. Ilie, W.F. Smyth.
 *   A simple algorithm for computing the Lempel-Ziv
 *   factorization, in: J.A. Storer, M.W. Marcellini
 *   (Eds.), 18th Data Compression Conference, IEEE
 *   Computer Society, Los Alamitos, CA, 2008,
 *   pp. 482-488.
 */
size_t *computeLpf(Esa *esa, int64_t **prevOccP) {
  size_t n = esa->n;
  saidx_t *sa = esa->sa;
  int64_t *lcp = esa->lcp;

  size_t *lpf = (size_t *)emalloc((n + 1) * sizeof(size_t));
  int64_t *prevOcc = (int64_t *)emalloc(n * sizeof(int64_t));

  lpf[n] = 0;
  kvec_t(Elem) s;
  kv_init(s);
  kv_push(Elem, s, ((Elem){0, 0}));

  for (size_t i = 1; i <= n; i++) {
    int64_t currLcp = i == n ? 0 : lcp[i];
    int64_t sai = i == n ? -1 : sa[i];
    while (!kv_empty(s) && sai < sa[kv_top(s).i]) {
      Elem v = kv_pop(s);
      lpf[sa[v.i]] = MAX(v.lcp, currLcp);
      currLcp = MIN(v.lcp, currLcp);
      // fill prevOcc
      if (lpf[sa[v.i]] == 0)
        prevOcc[sa[v.i]] = -1;
      else if (v.lcp > currLcp)
        prevOcc[sa[v.i]] = sa[kv_top(s).i];
      else
        prevOcc[sa[v.i]] = sai;
    }
    if (i < n)
      kv_push(Elem, s, ((Elem){i, currLcp}));
  }
  kv_destroy(s);

  *prevOccP = prevOcc;
  return lpf;
}

// alternative prevOcc calculation (from same paper)
size_t *computeLpf2(Esa *esa, int64_t **prevOccP) {
  size_t n = esa->n;
  saidx_t *sa = esa->sa;
  int64_t *lprev = emalloc(n * sizeof(int64_t));
  int64_t *lnext = emalloc(n * sizeof(int64_t));
  int64_t *prevl = emalloc(n * sizeof(int64_t));
  int64_t *prevr = emalloc(n * sizeof(int64_t));
  size_t *lpf = ecalloc(n, sizeof(int64_t));
  int64_t *lpfl = ecalloc(n, sizeof(int64_t));
  int64_t *lpfr = ecalloc(n, sizeof(int64_t));
  int64_t *prevOcc = emalloc(n * sizeof(int64_t));
  for (size_t i = 0; i < esa->n; i++) {
    lprev[sa[i]] = i == 0 ? -1 : sa[i - 1];
    lnext[sa[i]] = i == esa->n - 1 ? -1 : sa[i + 1];
  }
  for (int64_t sai = esa->n - 1; sai >= 0; sai--) {
    prevl[sai] = lprev[sai];
    prevr[sai] = lnext[sai];
    if (lprev[sai] != -1)
      lnext[lprev[sai]] = lnext[sai];
    if (lnext[sai] != -1)
      lprev[lnext[sai]] = lprev[sai];
  }
  free(lprev);
  free(lnext);

  prevOcc[0] = -1;
  for (size_t i = 1; i < esa->n; i++) {
    size_t j = MAX(lpfl[i - 1] - 1, 0);
    size_t k = MAX(lpfr[i - 1] - 1, 0);
    if (prevl[i] == -1)
      lpfl[i] = 0;
    else {
      while (esa->str[i + j] == esa->str[prevl[i] + j])
        j++;
      lpfl[i] = j;
    }
    if (prevr[i] == -1)
      lpfr[i] = 0;
    else {
      while (esa->str[i + k] == esa->str[prevr[i] + k])
        k++;
      lpfr[i] = k;
    }
    lpf[i] = MAX(lpfl[i], lpfr[i]);
    if (lpf[i] == 0)
      prevOcc[i] = -1;
    else if (lpfl[i] > lpfr[i])
      prevOcc[i] = prevl[i];
    else
      prevOcc[i] = prevr[i];
  }

  /* for (size_t i = 0; i < esa->n; i++) */
  /*   printf("i:%zu sai:%ld prevl:%ld prevr:%ld lpf:%ld po:%ld\n", */
  /*       i, sa[i], prevl[sa[i]], prevr[sa[i]], lpf[sa[i]],prevOcc[sa[i]]); */

  free(prevl);
  free(prevr);
  free(lpfl);
  free(lpfr);
  *prevOccP = prevOcc;
  return lpf;
}

Fact *computeLZFact(Esa *esa, bool alternative) {
  size_t n = esa->n;
  int64_t *prevOcc;
  size_t *lpf = alternative ? computeLpf2(esa, &prevOcc) : computeLpf(esa, &prevOcc);

  Fact *lzf = (Fact *)emalloc(sizeof(Fact));
  lzf->fact = (size_t *)emalloc(n * sizeof(size_t));
  lzf->fact[0] = 0;
  long i = 0;
  while (lzf->fact[i] < n) {
    lzf->fact[i + 1] = lzf->fact[i] + MAX(1, lpf[lzf->fact[i]]);
    i++;
  }
  lzf->lpf = lpf;
  lzf->prevOcc = prevOcc;
  lzf->n = i;

  lzf->str = esa->str;
  lzf->strLen = esa->n;
  return lzf;
}
