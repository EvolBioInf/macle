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

List **getPeriodicityLists(bool runsOnly, Fact *lzf, size_t **pnum);
Periodicity *collectPeriodicities(List **pl, size_t seqLen, size_t pnum);
void freePeriodicityLists(List **pl, size_t seqLen, size_t *pnum);

Periodicity *getPeriodicities2(Esa *esa, size_t *pnum);
