#pragma once
#include <vector>
#include <list>

#include "periodicity.h"

uint64_t maxFacts(uint64_t g, uint64_t n);

void mlComplexity(size_t w, size_t k, std::vector<double> &y, Fact &mlf, double gc);
void runComplexity(size_t w, size_t k, std::vector<double> &y, size_t n,
                   std::vector<std::list<Periodicity>> &ls);
