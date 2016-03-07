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
size_t *computeLpf(Esa *esa) {
  size_t n = esa->n;
  esa->sa = erealloc(esa->sa, (n + 1) * sizeof(uint64_t));
  esa->lcp = erealloc(esa->lcp, (n + 1) * sizeof(int64_t));
  size_t *lpf = (size_t *)emalloc((n + 1) * sizeof(size_t));

  uint64_t *sa = esa->sa;
  int64_t *lcp = esa->lcp;
  sa[n] = -1;
  lcp[0] = 0;
  lcp[n] = 0;
  lpf[n] = 0;
  Stack *s = newStack(1);
  stackPush(s, 0);

  for (size_t i = 1; i <= n; i++) {
    while (!stackEmpty(s) &&
           (sa[i] < sa[(size_t)stackTop(s)] ||
            (sa[i] > sa[(size_t)stackTop(s)] && lcp[i] <= lcp[(size_t)stackTop(s)]))) {
      if (sa[i] < sa[(size_t)stackTop(s)]) {
        lpf[sa[(size_t)stackTop(s)]] = MAX(lcp[(size_t)stackTop(s)], lcp[i]);
        lcp[i] = MIN(lcp[(size_t)stackTop(s)], lcp[i]);
      } else
        lpf[sa[(size_t)stackTop(s)]] = lcp[(size_t)stackTop(s)];
      stackPop(s);
    }
    if (i < n)
      stackPush(s, (void *)i);
  }
  freeStack(s);

  return lpf;
}

Fact *computeLZFact(Esa *esa) {
  size_t n = esa->n;
  size_t *lpf = computeLpf(esa);

  Fact *lzf = (Fact *)emalloc(sizeof(Fact));
  lzf->fact = (size_t *)emalloc(n * sizeof(size_t));
  lzf->fact[0] = 0;
  long i = 0;
  while (lzf->fact[i] < n) {
    lzf->fact[i + 1] = lzf->fact[i] + MAX(1, lpf[lzf->fact[i]]);
    i++;
  }
  lzf->lpf = lpf;
  lzf->n = i;

  lzf->str = esa->str;
  lzf->strLen = esa->n;
  return lzf;
}
