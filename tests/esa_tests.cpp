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

char const *test_getEsa() {
  FastaFile ff("Data/hotspotExample2.fasta");

  char *s = ff.seqs[0].seq;
  size_t n = ff.seqs[0].len;
  Esa esa(s, n + 1); // calculate esa, including $

  mu_assert(esa.str == s, "ESA does not point to original sequence");
  mu_assert(esa.n == n + 1, "ESA size not correct");
  mu_assert(esa.str[esa.sa[0]] == '$', "first ESA entry not $");
  mu_assert(esa.sa[0] == (int64_t)n, "wrong SA index");
  mu_assert(esa.isa[esa.sa[0]] == 0, "isa incorrect");
  mu_assert(esa.isa[esa.sa[n]] == (int64_t)n, "isa incorrect");
  mu_assert(esa.lcp[0] == -1, "first LCP not -1");
  mu_assert(esa.lcp[esa.n] == -1, "last LCP not -1");

  return NULL;
}

char const *test_revEsaRnd() {
  /* char *s = "AACCGGTTGGTT$"; // from Ohlebusch book */
  size_t n = 100;
  string str = randSeq(n++);
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
      /* printf("%zu %zu -> %zu %zu\n", i,j, n-i-2, n-j-2); */
      int64_t obs = resa.getLcp(rmq, n - i - 1, n - j - 1);
      if (exp != obs)
        printf("%zu %zu -> %zu %zu\n", i, j, n - i - 1, n - j - 1);
      mu_assert_eq(exp, obs, "lcs does not match");
    }

  return NULL;
}

char const *all_tests() {
  srand(time(NULL));
  mu_suite_start();
  mu_run_test(test_getEsa);
  mu_run_test(test_revEsaRnd);
  return NULL;
}

RUN_TESTS(all_tests)
