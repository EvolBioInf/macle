#pragma once
#include <vector>
#include <list>
#include <utility>

#include "periodicity.h"

size_t numEntries(size_t n, size_t w, size_t k);
std::vector<size_t> calcNAWindows(size_t n, size_t w, size_t k,
                                  std::vector<std::pair<size_t, size_t>> const &badiv);

void mlComplexity(size_t n, size_t w, size_t k, std::vector<double> &y,
                  std::vector<size_t> const &fact, double gc,
                  std::vector<size_t> const &badw = std::vector<size_t>(0));
void runComplexity(size_t n, size_t w, size_t k, std::vector<double> &y,
                   std::vector<std::list<Periodicity>> const &ls, double gc,
                   std::vector<size_t> const &badw = std::vector<size_t>(0), bool calcAvg=false);

double calcExpRunComplexity(size_t len, double gc, size_t reps);
