#pragma once

#include "esa.h"

typedef struct interval {
  int64_t lcp;
  size_t lb;
  size_t rb;

  struct interval **childList;
  size_t numChildren;
} Interval;

Interval *newInterval(int lcp, int lb, int rb);
void addChild(Interval *parent, Interval *child);
void printInterval(Interval *in);
void printLcpTree(Interval *in);

bool getSubInterval(Interval *ret, Esa *esa, Interval iv, char c);
Interval getInterval(Esa *esa, char *query, size_t n);

Interval *getLcpTree(Esa *esa);
void freeLcpTree(Interval *iv);
int64_t getLcpWithTree(Esa *esa, Interval *tree, size_t i, size_t j);
