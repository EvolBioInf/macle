/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include "matchlength.h"
#include "shulen.h"
#include <algorithm>

// calculate number of factor start positions up to some given prefix length
// . allows to easily get number of factor start positions in any interval
uint64_t *computeFactPrefixSum(Fact &f) {
  size_t n = f.strLen;
  /* compute observed number of match factors for every prefix */
  uint64_t *ps = new uint64_t[n];
  size_t nextfact = 1;
  ps[0] = 1;
  for (size_t i = 1; i < n; i++) {
    ps[i] = ps[i - 1];
    if (nextfact < f.n && i == f.fact[nextfact]) {
      ps[i]++;
      nextfact++;
    }
  }
  return ps;
}

Fact computeMLFact(Esa &esa) {
  Fact mlf;
  mlf.prevOcc = nullptr; // no prevOcc array here
  mlf.str = esa.str;
  mlf.strLen = esa.n;

  /* construct and fill array of match lengths */
  uint64_t *ml = new uint64_t[esa.n];
  for (size_t i = 0; i < esa.n; i++) {
    ml[esa.sa[i]] = std::max((int64_t)1, std::max(esa.lcp[i], esa.lcp[i + 1]));
  }

  /* compute observed number of match factors, store their positions */
  mlf.fact = new size_t[esa.n];
  mlf.n = 0;
  size_t i = 0;
  while (i < esa.n) {
    mlf.fact[mlf.n++] = i;
    i += ml[i];
  }
  delete[] ml;

  // use lpf field to store prefix sums for ML factors
  mlf.lpf = computeFactPrefixSum(mlf);

  return mlf;
}
