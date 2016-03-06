#pragma once
#include "esa.h"

typedef struct Periodicity {
  size_t b; //beginning index
  size_t e; //end index (inclusive)
  size_t l; //period length
} Periodicity;

void printPeriodicity(Periodicity *p);

Periodicity *getPeriodicities(Fact *lzf, size_t *plen);
Periodicity *getPeriodicities2(char *str, size_t n, size_t *plen);

