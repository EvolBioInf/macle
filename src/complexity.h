#pragma once
#include <vector>
#include <list>
#include <utility>

#include "index.h"
#include "periodicity.h"

size_t numEntries(size_t n, size_t w, size_t k);

void mlComplexity(size_t offset, size_t n, size_t w, size_t k, std::vector<double> &y, ComplexityData const &dat);
void runComplexity(size_t offset, size_t n, size_t w, size_t k, std::vector<double> &y, ComplexityData const &dat, bool calcAvg=true);

double calcAvgRunComplexity(size_t len, double gc, size_t reps);
double estimateAvgRunComplexity(double gc);
