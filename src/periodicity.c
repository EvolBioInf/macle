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

// TODO: maybe efficient internal lcp, and lcs (longest common suffix) (Ohlebusch!)

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
        ps[*plen].b = i - l - L - 1; // again 0-indexed
        ps[*plen].e = i - 1 + R - 1;
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

// helper to avoid duplication
static inline bool addPeriodicity(bool runsOnly, List **Lt1, size_t b, size_t e,
                                  size_t l) {
  if (Lt1[b] && runsOnly) { // if we want only runs,
    List *last = listLast(Lt1[b]);
    if (((Periodicity *)last->value)->e != e) { // add only if prev. is not same interval
      listAppend(&last, newPeriodicity(b, e, l));
      return true;
    }
  } else { // otherwise, add in any case
    listAppend(&Lt1[b], newPeriodicity(b, e, l));
    return true;
  }
  return false;
}

// Algorithm 5.17 - calculate type1 periodicities
List **calcType1Periodicities(bool runsOnly, Fact *lzf, size_t *pnum) {
  (*pnum) = 0;
  // as required, array of lists of max. per. of type 1
  // indexed by start position, each sorted by end position
  List **Lt1 = ecalloc(lzf->strLen + 1, sizeof(List *));
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
        size_t b = bj - l - L - 1;
        size_t e = ejm1 + R - 1;
        if (addPeriodicity(runsOnly, Lt1, b, e, l))
          (*pnum)++;
      }
    }

    for (size_t l = 1; l <= sjLen; l++) {
      size_t L = lcs2(lzf->str, ejm1, ejm1 + l);
      size_t R = lcp2(lzf->str, lzf->strLen, bj, bj + l);
      if (L + R >= l && bj + l - 1 + R <= ej && L < l) {
        size_t b = bj - L - 1;
        size_t e = ejm1 + l + R - 1;
        if (addPeriodicity(runsOnly, Lt1, b, e, l))
          (*pnum)++;
      }
    }
  }

  return Lt1;
}

// prevOcc conforming to algorithm: 1-indexed and current value instead of -1
static inline int64_t getPrevOcc(Fact *lzf, size_t i) {
  int64_t tmp = lzf->prevOcc[lzf->fact[i]];
  return (tmp == -1 ? (int64_t)lzf->fact[i] : tmp) + 1;
}

// add given periodicity into corresponding same-beginning list
// so that list remains sorted by end positions. This hack is necessary!
// Ohlebusch/Algorithm 5.18 says "prepend". This is WRONG!!! Was a fun night
// debugging -.-'
void listInsert(List **l, Periodicity *p) {
  List *curr = *l;
  if (!curr || ((Periodicity *)curr->value)->e >= p->e) {
    listPrepend(l, p);
    return;
  }
  while (curr->next && ((Periodicity *)curr->next->value)->e < p->e)
    curr = curr->next;
  List *item = newList();
  item->value = p;
  item->next = curr->next;
  curr->next = item;
}

// Algorithm 5.18 - calculate type 2 periodicities (proper substrings of LZ-factors)
void calcType2Periodicities(List **Lt1, Fact *lzf, size_t *pnum) {
  for (size_t j = 1; j < lzf->n; j++) {
    size_t bj = factStart(lzf, j); // b_j
    size_t ej = factEnd(lzf, j);   // e_j
    if (factLen(lzf, j) >= 4) {
      int64_t dj = bj - getPrevOcc(lzf, j);
      for (size_t i = bj; i < ej - 1; i++) {
        for (eachListItem(curr, Lt1[i - dj])) {
          Periodicity *p = (Periodicity *)(curr->value);
          if (p->e + dj >= ej - 1)
            break;
          Periodicity *newp = newPeriodicity(i, p->e + dj, p->l);
          /* listPrepend(&Lt1[i], newp); */ // XXX: doesn't work with prepend!!!
          listInsert(&Lt1[i], newp); // need to sort into list correctly!!!!!!!!!!!!!!!
          (*pnum)++;
        }
      }
    }
  }
}

// max. periodicities in O(n) using LZ-Factors, returns array of lists
List **getPeriodicityLists(bool runsOnly, Fact *lzf, size_t *pnum) {
  List **Lt1 = calcType1Periodicities(runsOnly, lzf, pnum);
  calcType2Periodicities(Lt1, lzf, pnum);
  return Lt1;
}

Periodicity *collectPeriodicities(List **pl, size_t seqLen, size_t pnum) {
  Periodicity *ps = emalloc(pnum * sizeof(Periodicity));
  size_t ind = 0;
  for (size_t i = 0; i < seqLen; i++) {
    List *curr = pl[i];
    List *tmp;
    if (curr) {
      do {
        Periodicity *p = (Periodicity *)(curr->value);
        ps[ind++] = *p;
        free(p);
        tmp = curr;
        curr = curr->next;
        free(tmp);
      } while (curr);
    }
  }
  free(pl);
  return ps;
}

void freePeriodicityLists(List **pl, size_t seqLen) {
  for (size_t i = 0; i < seqLen; i++) {
    List *curr = pl[i];
    List *tmp;
    if (curr) {
      do {
        free((Periodicity *)(curr->value));
        tmp = curr;
        curr = curr->next;
        free(tmp);
      } while (curr);
    }
  }
  free(pl);
}

Periodicity *getPeriodicities(bool runsOnly, Fact *lzf, size_t *pnum) {
  List **pl = getPeriodicityLists(runsOnly, lzf, pnum);
  return collectPeriodicities(pl, lzf->strLen, *pnum);
}
