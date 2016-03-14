#pragma once
#include "esa.h"
#include "factors.h"
#include "list.h"

typedef struct Periodicity {
  size_t b; // beginning index
  size_t e; // end index (inclusive)
  size_t l; // period length
} Periodicity;

void printPeriodicity(Periodicity *p);

// these only exported for test cases
int64_t lcp(Esa *esa, int64_t *lcptab, size_t i, size_t j);
int64_t lcs(Esa *resa, int64_t *rlcptab, size_t i, size_t j);
size_t lcs2(char *str, int64_t i, int64_t j);
size_t lcp2(char *str, size_t n, size_t i, size_t j);

Periodicity *getPeriodicities(bool runsOnly, Fact *lzf, Esa *esa, size_t *pnum);

List **getPeriodicityLists(bool runsOnly, Fact *lzf, Esa *esa, size_t *pnum);
Periodicity *collectPeriodicities(List **pl, size_t seqLen, size_t pnum);
void freePeriodicityLists(List **pl, size_t seqLen);

Periodicity *getPeriodicities2(Esa *esa, size_t *pnum);

/* number of periodicities corresponding to a run */
static inline size_t persFromRun(Periodicity *p) {
  return (p->e - p->b + 1) / (2 * p->l);
}
