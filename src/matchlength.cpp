/***** complexity.c *******************************
 * Description:
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Jul 15 10:49:56 2015
 **************************************************/
#include "matchlength.h"
#include "shulen.h"
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;

void Fact::print() const {
  stringstream ss;
  size_t n = fact.size();
  for (size_t i = 0; i < n; i++) {
    size_t start = this->fact[i];
    size_t end = i < n - 1 ? this->fact[i + 1] : this->strLen;
    ss << string(this->str + start, end-start) << (i < n - 1 ? "." : "");
  }
  string s = ss.str();
  for (size_t i = 0; i < s.size(); i+=80)
    cout << s.substr(i, min(80UL, s.size()-i)) << endl;
}

// length of factor
size_t Fact::factLen(size_t i) const {
  Fact const &f = *this;
  if (i == 0)
    return f.fact[1];
  if (i == f.fact.size() - 1)
    return f.strLen - f.fact[f.fact.size() - 1];
  return f.fact[i + 1] - f.fact[i];
}

//input: esa for both strands (seq$revcompseq$)
void computeMLFact(Fact &mlf, Esa const &esa) {
  mlf.fact.resize(0);
  mlf.str = esa.str;
  mlf.strLen = esa.n/2; //single strand length

  /* construct and fill array of match lengths */
  vector<uint64_t> ml(esa.n);
  for (size_t i = 0; i < esa.n; i++) {
    ml[esa.sa[i]] = max(1UL, (size_t)max(esa.lcp[i], esa.lcp[i + 1]));
  }

  /* compute observed number of match factors, store their positions */
  vector<size_t> factmp;
  size_t i = 0;
  while (i < mlf.strLen) {
    factmp.push_back(i);
    i += ml[i];
  }

  mlf.fact.resize(factmp.size());
  for (i=0; i<factmp.size(); i++)
    mlf.fact[i] = factmp[i];
#ifdef USE_SDSL
    sdsl::util::bit_compress(mlf.fact);
#endif
}
