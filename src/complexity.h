#pragma once
#include <vector>
#include <list>

#include "periodicity.h"

void mlComplexity(size_t w, size_t k, std::vector<double> &y, size_t n,
                  std::vector<size_t> const &fact, double gc);
void runComplexity(size_t w, size_t k, std::vector<double> &y, size_t n,
                   std::vector<std::list<Periodicity>> const &ls);
