#include "minunit.h"
#include <string.h>
#include <time.h>

#include "sequenceData.h"
#include "stringUtil.h"
#include "esa.h"

// naive: length of lcs of prefixes 1..i and 1..j of given seq (input 1-indexed)
size_t lcsNaive(char *str, int64_t i, int64_t j) {
  int64_t k = 0;
  while (((int64_t)MIN(i, j)) - k >= 0 && str[i - k] == str[j - k])
    k++;
  return k;
}

char *test_getEsa() {
  Sequence *seq = readFastaFromFile("Data/hotspotExample2.fasta");

  char *s = seqStr(seq, 0);
  size_t n = seqLen(seq, 0);
  Esa *esa = getEsa(s, n + 1); // calculate esa, including $

  mu_assert(esa->str == s, "ESA does not point to original sequence");
  mu_assert(esa->n == n + 1, "ESA size not correct");
  mu_assert(esa->str[esa->sa[0]] == '$', "first ESA entry not $");
  mu_assert(esa->sa[0] == (int64_t)n, "wrong SA index");
  mu_assert(esa->isa[esa->sa[0]] == 0, "isa incorrect");
  mu_assert(esa->isa[esa->sa[n]] == (int64_t)n, "isa incorrect");
  mu_assert(esa->lcp[0] == -1, "first LCP not -1");
  mu_assert(esa->lcp[esa->n] == -1, "last LCP not -1");

  freeEsa(esa);
  freeSequence(seq);
  return NULL;
}

char *test_revEsaRnd() {
  /* char *s = "AACCGGTTGGTT$"; // from Ohlebusch book */
  size_t n = 100;
  char *s = randSeq(n++);
  Esa *esa = getEsa(s, n); // calculate esa, including $

  char *srev = strdup2(esa->str);
  reverse(srev, n);
  Esa *resa = getEsa(srev, n);
  int64_t *revlcptab = precomputeLcp(resa);

  /* printEsa(esa); */
  /* printEsa(resa); */

  for (size_t i = 0; i < n - 1; i++)
    for (size_t j = i + 1; j < n; j++) {
      int64_t exp = lcsNaive(s, i, j);
      /* printf("%zu %zu -> %zu %zu\n", i,j, n-i-2, n-j-2); */
      int64_t obs = getLcp(resa, revlcptab, n - i - 1, n - j - 1);
      if (exp != obs)
        printf("%zu %zu -> %zu %zu\n", i, j, n - i - 1, n - j - 1);
      mu_assert_eq(exp, obs, "lcs does not match");
    }

  freeEsa(resa);
  freeEsa(esa);
  free(revlcptab);
  free(srev);
  free(s);
  return NULL;
}

char *all_tests() {
  srand(time(NULL));
  mu_suite_start();
  mu_run_test(test_getEsa);
  mu_run_test(test_revEsaRnd);
  return NULL;
}

RUN_TESTS(all_tests)
