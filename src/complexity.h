#pragma once

uint64_t factorsFromTo(Fact *f, int64_t l, int64_t r);

void mlComplexity(size_t w, size_t k, double *y, Fact *mlf, double gc);
void runComplexity(size_t w, size_t k, double *y, size_t n, List **ls, size_t *lens);
