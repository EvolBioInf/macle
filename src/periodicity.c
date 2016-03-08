#include "prelude.h"
#include "eprintf.h"

#include "lempelziv.h"
#include "periodicity.h"
#include "esa.h"
#include "interval.h"
#include "list.h"

Periodicity *newPeriodicity(size_t b, size_t e, size_t l) {
  Periodicity *p = emalloc(sizeof(Periodicity));
  p->b = b;
  p->e = e;
  p->l = l;
  return p;
}

void printPeriodicity(Periodicity *p) { printf("(%zu,%zu,%zu)\n", p->b, p->e, p->l); }

// TODO: efficient internal lcp, and lcs (longest common suffix) (Ohlebusch!)

// naive: length of lcp of suffixes i and j of given seq (1-indexed)
size_t lcp2(char *str, size_t n, size_t i, size_t j) {
  i--;
  j--;
  uint64_t k = 0;
  while (MAX(i, j) + k < n && str[k + i] == str[k + j])
    k++;
  return k;
}

// naive: length of lcs of prefixes 1..i and 1..j of given seq (input 1-indexed)
size_t lcs2(char *str, int64_t i, int64_t j) {
  i--;
  j--;
  int64_t k = 0;
  while (((int64_t)MIN(i, j)) - k >= 0 && str[i - k] == str[j - k])
    k++;
  return k;
}

// max. periodicities in O(nlog(n)) (Algorithm 5.16), mainly for reference / comparison
Periodicity *getPeriodicities2(Esa *esa, size_t *plen) {
  size_t n = esa->n;
  char *str = esa->str;
  Periodicity *ps = emalloc(((size_t)2.05 * n) * sizeof(Periodicity));
  *plen = 0;
  for (size_t l = 1; l <= n / 2; l++) {
    size_t i = 2 * l + 1;
    while (i <= n + 1) {
      size_t L = lcs2(str, i - 1 - l, i - 1);
      size_t R = lcp2(str, n, i - l, i);
      if (L + R >= l && R < l) {
        ps[*plen].b = i - l - L;
        ps[*plen].e = i - 1 + R;
        ps[*plen].l = l;
        (*plen)++;
      }
      i += l;
    }
  }
  return ps;
}

// start position (1-indexed)
static inline size_t factStart(Fact *f, size_t i) { return f->fact[i] + 1; }

// end position (1-indexed)
static inline size_t factEnd(Fact *f, size_t i) {
  return (i < f->n - 1 ? f->fact[i + 1] : f->strLen);
}

// Algorithm 5.17 - calculate type1 periodicities
List **calcType1Periodicities(bool runsOnly, Fact *lzf) {
  // as required, array of lists of max. per. of type 1
  // indexed by start position, each sorted by end position
  List **Lt1 = ecalloc(lzf->strLen, sizeof(List *));
  for (size_t j = 1; j < lzf->n; j++) {
    size_t sjLen = factLen(lzf, j);       // |s_j|
    size_t sjm1Len = factLen(lzf, j - 1); // |s_(j-1)|
    size_t bj = factStart(lzf, j);        // b_j
    size_t bjm1 = factStart(lzf, j - 1);  // b_(j-1)
    size_t ej = factEnd(lzf, j);          // e_j
    size_t ejm1 = factEnd(lzf, j - 1);    // e_(j-1)

    size_t max = MIN(sjLen + sjm1Len - 1, ejm1);
    for (size_t l = 1; l <= max; l++) {
      size_t L = lcs2(lzf->str, ejm1 - l, ejm1);
      size_t R = lcp2(lzf->str, lzf->strLen, bj - l, bj);
      if (L + R >= l && (R >= 1 || bj - l - L > bjm1)) {
        Periodicity *p = newPeriodicity(bj - l - L, ejm1 + R, l);
        listAppend(&Lt1[p->b], p);
      }
    }

    for (size_t l = 1; l <= sjLen; l++) {
      size_t L = lcs2(lzf->str, ejm1, ejm1 + l);
      size_t R = lcp2(lzf->str, lzf->strLen, bj, bj + l);
      if (L + R >= l && bj + l - 1 + R <= ej && L < l) {
        Periodicity *p = newPeriodicity(bj - L, ejm1 + l + R, l);
        listAppend(&Lt1[p->b], p);
      }
    }
  }

  if (runsOnly) { // remove out non-minimal periodicities -> keep "runs"
    for (size_t i = 0; i < lzf->strLen; i++) {
      Periodicity *lastp = NULL;
      List *last = NULL;
      List *curr = Lt1[i];
      if (curr)
        do {
          Periodicity *currp = (Periodicity *)curr->value;
          if (lastp && lastp->b == currp->b && lastp->e == currp->e) {
            last->next = curr->next;
            free(currp);
            free(curr);
            curr = last->next;
          } else {
            last = curr;
            curr = curr->next;
            lastp = currp;
          }
        } while (curr);
    }
  }
  return Lt1;
}

// calculate PrevOcc_j (previous occurence of LZ-factor s_j in text, 1-indexed)
// if a factor is a first occurence, returns -1
// TODO: use fact that multiple s_j can be the same -> organize as trie,
// to reuse lcp-intervals for calculation! MUCH TOO SLOW
size_t *calcPrevOcc(Fact *lzf, Esa *esa) {
  size_t *prevOcc = emalloc(lzf->n * sizeof(size_t));

  for (size_t i = 0; i < lzf->n; i++) {
    size_t len = factLen(lzf, i);
    /* printf("getInterval for %zu/%zu (len: %zu)\n", i+1, lzf->n, len); */
    Interval iv = getInterval(esa, lzf->str + lzf->fact[i], len);
    int64_t best = -1;
    for (size_t j = iv.lb; j <= (size_t)iv.rb; j++) {
      size_t curr = esa->sa[j];
      if (curr < lzf->fact[i] && (int64_t)curr > best)
        best = curr;
    }
    prevOcc[i] = best + 1; // 1-indexed!
  }

  return prevOcc;
}

// Algorithm 5.18 - calculate type 2 periodicities (proper substrings of LZ-factors)
void calcType2Periodicities(List **Lt1, Fact *lzf, Esa *esa) {
  size_t *prevOcc = calcPrevOcc(lzf, esa);
  for (size_t j = 1; j < lzf->n; j++) {
    size_t bj = factStart(lzf, j); // b_j
    size_t ej = factEnd(lzf, j);   // e_j
    if (factLen(lzf, j) >= 4) {
      size_t dj = bj - prevOcc[j];
      for (size_t i = bj + 1; i <= ej - 1; i++) {
        List *curr = Lt1[i - dj];
        if (curr)
          do {
            Periodicity *p = (Periodicity *)(curr->value);
            if (p->e + dj >= ej)
              break;
            Periodicity *newp = newPeriodicity(i, p->e + dj, p->l);
            listPrepend(&Lt1[i], newp);
            curr = curr->next;
          } while (curr);
      }
    }
  }
  free(prevOcc);
}

// max. periodicities in O(n) using LZ-Factors
Periodicity *getPeriodicities(bool runsOnly, Fact *lzf, Esa *esa, size_t *plen) {
  List **Lt1 = calcType1Periodicities(runsOnly, lzf);
  calcType2Periodicities(Lt1, lzf, esa);

  // collect results from lists and cleanup
  // Theorem 5.3.[39,43] (upper bound for num. max. per. / runs for a string)
  size_t num = (size_t)(lzf->strLen * (runsOnly ? 1.03 : 2.05));
  Periodicity *ps = emalloc(num * sizeof(Periodicity));
  *plen = 0;
  for (size_t i = 0; i < lzf->strLen; i++) {
    List *curr = Lt1[i];
    List *tmp;
    if (curr) {
      do {
        Periodicity *p = (Periodicity *)(curr->value);
        ps[(*plen)++] = *p;
        free(p);
        tmp = curr;
        curr = curr->next;
        free(tmp);
      } while (curr);
    }
  }
  free(Lt1);
  return ps;
}
