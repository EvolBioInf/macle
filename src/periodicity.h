#pragma once
#include "esa.h"
#include "factors.h"

#include <vector>
#include <list>

struct Periodicity {
  Periodicity(size_t b, size_t e, size_t l);
  size_t b; // beginning index
  size_t e; // end index (inclusive)
  size_t l; // period length
};

void printPeriodicity(Periodicity &p);

// these only exported for test cases
int64_t lcp(Esa &esa, int64_t *lcptab, size_t i, size_t j);
int64_t lcs(Esa &resa, int64_t *rlcptab, size_t i, size_t j);
size_t lcs2(char const *str, int64_t i, int64_t j);
size_t lcp2(char const *str, size_t n, size_t i, size_t j);

std::vector<Periodicity> getPeriodicities(bool runsOnly, Fact &lzf, Esa &esa,
                                          size_t &pnum);

std::vector<std::list<Periodicity>> getPeriodicityLists(bool runsOnly, Fact &lzf,
                                                        Esa &esa, size_t &pnum);
std::vector<Periodicity> collectPeriodicities(std::vector<std::list<Periodicity>> &lst);

std::vector<Periodicity> getPeriodicities2(Esa &esa);

/* length (from start to end) of run/periodicity */
static inline size_t perLen(Periodicity &p) { return p.e - p.b + 1; }

/* number of periodicities corresponding to a run */
static inline size_t persFromRun(Periodicity &p) { return perLen(p) / (2 * p.l); }
