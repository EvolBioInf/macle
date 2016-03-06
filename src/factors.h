#pragma once
#include "esa.h"

/* a factorization of a string (LempelZiv or MatchLength) */
typedef struct Fact{
  size_t *fact;   /* positions of factors */
  size_t *lpf;    /* lpf in case of Lempel-Ziv, otherwise unset */
  size_t n;       /* number of factors */

  char *str;      /* string */
  size_t strLen;  /* string length */

  double cObs, cMax, cMin, cNor;
} Fact;

void freeFact(Fact *f);
void printFact(Fact *f);
