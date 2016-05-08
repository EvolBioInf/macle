#pragma once
#include <vector>
#include <list>

#include "periodicity.h"

size_t numEntries(size_t n, size_t w, size_t k);

void mlComplexity(size_t n, size_t w, size_t k, std::vector<double> &y,
                  std::vector<size_t> const &fact, double gc);
void runComplexity(size_t n, size_t w, size_t k, std::vector<double> &y,
                   std::vector<std::list<Periodicity>> const &ls);

// runs that are smaller that this are thrown away -> seen as random noise
#define FILTER 10
