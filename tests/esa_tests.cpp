#include "minunit.h"
#include "util.h"
#include <cstring>
#include <ctime>
#include <algorithm>
#include <string>
using namespace std;

#include "fastafile.h"
#include "esa.h"

// naive: length of lcs of prefixes 1..i and 1..j of given seq (input 1-indexed)
size_t lcsNaive(char const *str, int64_t i, int64_t j) {
  int64_t k = 0;
  while (((int64_t)min(i, j)) - k >= 0 && str[i - k] == str[j - k])
    k++;
  return k;
}

// some basic tests
void test_getEsa() {
  FastaFile ff("Data/hotspotExample2.fasta");
  ff.seqs[0].seq += "$";

  char const *s = ff.seqs[0].seq.c_str();
  size_t n = ff.seqs[0].seq.size();
  Esa esa(s, n); // calculate esa, including $

  mu_assert(esa.str == s, "ESA does not point to original sequence");
  mu_assert(esa.n == n, "ESA size not correct");
  mu_assert(esa.str[esa.sa[0]] == '$', "first ESA entry not $");
  mu_assert(esa.sa[0] == (int64_t)(n - 1), "wrong SA index");
  mu_assert(esa.isa[esa.sa[0]] == 0, "isa incorrect");
  mu_assert(esa.isa[esa.sa[n - 1]] == (int64_t)(n - 1), "isa incorrect");
  mu_assert(esa.lcp[0] == -1, "first LCP not -1");
  mu_assert(esa.lcp[esa.n] == -1, "last LCP not -1");
}

// tests lcs value retrieval from reverse esa using getLcp vs naive
void test_revEsaRnd() {
  /* char *s = "AACCGGTTGGTT$"; // from Ohlebusch book */
  size_t n = 100;
  string str = randSeq(n++);
  str += "$";
  char const *s = str.c_str();
  Esa esa(s, n); // calculate esa, including $

  string strrev = str;
  reverse(strrev.begin(), strrev.end());
  char const *srev = strrev.c_str();

  Esa resa(srev, n);
  RMQ rmq = resa.precomputeLcp();

  /* printEsa(esa); */
  /* printEsa(resa); */

  for (size_t i = 0; i < n - 1; i++)
    for (size_t j = i + 1; j < n; j++) {
      int64_t exp = lcsNaive(s, i, j);
      int64_t obs = resa.getLcp(rmq, n - i - 1, n - j - 1);
      if (exp != obs)
        printf("%zu %zu -> %zu %zu\n", i, j, n - i - 1, n - j - 1);
      mu_assert_eq(exp, obs, "lcs does not match");
    }
}

// tests ESA reduction from seq$revcompseq$ -> seq$ without recalculation
void test_reduceEsa() {
  string str = randSeq(1000);
  string str2n = str + "$" + revComp(str) + "$";
  str += "$";
  Esa esaOne(str.c_str(), str.size());
  Esa esaBoth(str2n.c_str(), str2n.size());
  reduceEsa(esaBoth);
  esaBoth.str = str.c_str();

  mu_assert_eq(-1, esaBoth.lcp[esaBoth.sa.size()], "last LCP not -1!");
  for (size_t i = 0; i < str.size(); i++) {
    mu_assert_eq(esaOne.sa[i], esaBoth.sa[i], "SA[" << i << "] does not match");
    mu_assert_eq(esaOne.lcp[i], esaBoth.lcp[i], "LCP[" << i << "] does not match");
    mu_assert_eq(esaOne.isa[i], esaBoth.isa[i], "ISA[" << i << "] does not match");
  }
}

void all_tests() {
  srand(time(NULL));
  mu_run_test(test_getEsa);
  mu_run_test(test_revEsaRnd);
  mu_run_test(test_reduceEsa);
}
RUN_TESTS(all_tests)
