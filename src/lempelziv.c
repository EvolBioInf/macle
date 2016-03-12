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
#include "stack.h"
#include "eprintf.h"

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
  esa->sa = erealloc(esa->sa, (n + 1) * sizeof(int64_t));
  esa->lcp = erealloc(esa->lcp, (n + 1) * sizeof(int64_t));
  size_t *lpf = (size_t *)emalloc((n + 1) * sizeof(size_t));
  int64_t *prevOcc = (int64_t *)emalloc(n * sizeof(int64_t));

  int64_t *sa = esa->sa;
  int64_t *lcp = esa->lcp;
  sa[n] = -1;
  int64_t lcpn = lcp[n];
  lcp[n] = 0; // required for algorithm

  lpf[n] = 0;
  Stack *s = newStack(1);
  Stack *sLcp = newStack(1);
  stackPush(s, 0);
  stackPush(sLcp, 0);

  for (size_t i = 1; i <= n; i++) {
    int64_t currLcp = lcp[i];
    while (
        !stackEmpty(s) &&
        (sa[i] <
         sa[(size_t)stackTop(s)] /* ||  //TODO: clarify - this or-branch seems wrong?
            (sa[i] > sa[(size_t)stackTop(s)] && lcp[i] <= lcp[(size_t)stackTop(s)]) */)) {
      if (sa[i] < sa[(size_t)stackTop(s)]) {
        lpf[sa[(size_t)stackTop(s)]] = MAX((int64_t)stackTop(sLcp), currLcp);
        currLcp = MIN((int64_t)stackTop(sLcp), currLcp);
      } else
        lpf[sa[(size_t)stackTop(s)]] = (int64_t)stackTop(sLcp);
      int64_t v = (int64_t)stackPop(s);
      int64_t vLcp = (int64_t)stackPop(sLcp);
      // fill prevOcc
      if (lpf[sa[v]] == 0)
        prevOcc[sa[v]] = -1;
      else if (vLcp > currLcp)
        prevOcc[sa[v]] = sa[(size_t)stackTop(s)];
      else
        prevOcc[sa[v]] = sa[i];
    }
    if (i < n) {
      stackPush(s, (void *)i);
      stackPush(sLcp, (void *)currLcp);
    }
  }
  freeStack(s);
  freeStack(sLcp);
  lcp[n] = lcpn; // restore last value

  *prevOccP = prevOcc;
  return lpf;
}

// alternative prevOcc calculation (from same paper)
size_t *computeLpf2(Esa *esa, int64_t **prevOccP) {
  size_t n = esa->n;
  int64_t *sa = esa->sa;
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
