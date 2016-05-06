/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include "matchlength.h"
#include "shulen.h"
#include <algorithm>
#include <vector>
using namespace std;

void computeMLFact(Fact &mlf, Esa const &esa) {
  mlf.prevOcc.clear();
  mlf.fact.clear();
  mlf.lpf.clear();
  mlf.str = esa.str;
  mlf.strLen = esa.n;

  /* construct and fill array of match lengths */
  vector<uint64_t> ml(esa.n);
  for (size_t i = 0; i < esa.n; i++) {
    ml[esa.sa[i]] = std::max((int64_t)1, std::max(esa.lcp[i], esa.lcp[i + 1]));
  }

  /* compute observed number of match factors, store their positions */
  size_t i = 0;
  while (i < esa.n) {
    mlf.fact.push_back(i);
    i += ml[i];
  }
}
