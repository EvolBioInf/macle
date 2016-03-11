#pragma once
#include "esa.h"
#include "list.h"

typedef struct Periodicity {
  size_t b; // beginning index
  size_t e; // end index (inclusive)
  size_t l; // period length
} Periodicity;

void printPeriodicity(Periodicity *p);

Periodicity *getPeriodicities(bool runsOnly, Fact *lzf, size_t *pnum);

List **getPeriodicityLists(bool runsOnly, Fact *lzf, size_t *pnum);
Periodicity *collectPeriodicities(List **pl, size_t seqLen, size_t pnum);
void freePeriodicityLists(List **pl, size_t seqLen);

Periodicity *getPeriodicities2(Esa *esa, size_t *pnum);

/* number of periodicities corresponding to a run */
static inline size_t persFromRun(Periodicity *p) {
  return (p->e - p->b + 1) / (2 * p->l);
}
