#pragma once
#include "factors.h"
#include "esa.h"

Fact *computeMLFact(Esa *esa);
Fact *mlComplexity(Esa *esa, double gc);
uint64_t factorsFromTo(Fact *f, int64_t l, int64_t r);
