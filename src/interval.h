#pragma once

#include "esa.h"

typedef struct interval {
  int64_t lcp;
  int64_t lb;
  int64_t rb;

  struct interval **childList;
  size_t numChildren;
} Interval;

Interval *newInterval(int lcp, int lb, int rb);
void addChild(Interval *parent, Interval *child);
void printInterval(Interval *in);

bool getSubInterval(Interval *ret, Esa *esa, Interval iv, char c);
Interval getInterval(Esa *esa, char *query, size_t n);

// TODO: fix these
Interval *getLcpTree(Esa *esa);
int64_t getLcp(Interval *tree, size_t i, size_t j);
