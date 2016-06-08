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
int64_t lcp(const Esa &esa, const sdsl::rmq_succinct_sct<> &tab, size_t i, size_t j);
int64_t lcs(const Esa &resa, const sdsl::rmq_succinct_sct<> &rtab, size_t i, size_t j);
size_t lcs2(char const *str, size_t i, size_t j);
size_t lcp2(char const *str, size_t n, size_t i, size_t j);

std::vector<Periodicity> getPeriodicities(bool runsOnly, Fact const &lzf, Esa const &esa,
                                          size_t &pnum);

typedef std::vector<std::vector<Periodicity>> PerLists;

PerLists getPeriodicityLists(bool runsOnly, Fact const &lzf, Esa const &esa, size_t &pnum);
std::vector<Periodicity> collectPeriodicities(PerLists &lst);

std::vector<Periodicity> getPeriodicities2(Esa const &esa);

/* length (from start to end) of run/periodicity */
static inline size_t perLen(Periodicity const &p) { return p.e - p.b + 1; }

// number of periodicities corresponding to a run
/*
static inline size_t persFromRun(Periodicity const &p) { return perLen(p) / (2 * p.l); }
*/
