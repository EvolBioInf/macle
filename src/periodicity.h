#pragma once
#include "esa.h"

typedef struct Periodicity {
  size_t b; // beginning index
  size_t e; // end index (inclusive)
  size_t l; // period length
} Periodicity;

void printPeriodicity(Periodicity *p);

Periodicity *getPeriodicities(bool runsOnly, Fact *lzf, Esa *esa, size_t *plen);
Periodicity *getPeriodicities2(Esa *esa, size_t *plen);
