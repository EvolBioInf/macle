/***** complexity.h *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:55:13 2015
 **************************************************/
#pragma once
#include "factors.h"
#include "esa.h"

Fact *computeMLFact(Esa *esa);
Fact *mlComplexity(Esa *esa, double gc);
uint64_t factorsFromTo(Fact *f, int64_t l, int64_t r);
