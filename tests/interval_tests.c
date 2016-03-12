#include <time.h>
#include <string.h>
#include "minunit.h"

#include "stringUtil.h"
#include "esa.h"
#include "interval.h"

// hotspot paper example
static char *seq = "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"
                   "CACATATGCTAACTCTCAGTCTGTGTGTGCA$";

static char *seq2 = "AACCAACCAACCAA$"; // from Ohlebusch book

char *test_intervals() {
  Esa *esa = getEsa(seq, strlen(seq)); // calculate esa, including $

  Interval iv;
  // do some perfect matches
  iv = getInterval(esa, "$", 1);
  mu_assert(iv.lb == 0 && iv.rb == 0 && iv.lcp == 1, "wrong interval for $");

  iv = getInterval(esa, "AA", 2);
  mu_assert(iv.lb == 2 && iv.rb == 2 && iv.lcp == 2, "wrong interval for AA");

  iv = getInterval(esa, "ACGC", 4);
  mu_assert(iv.lb == 34 && iv.rb == 35 && iv.lcp == 4, "wrong interval for ACGC");

  iv = getInterval(esa, "CACA", 4);
  mu_assert(iv.lb == 41 && iv.rb == 71 && iv.lcp == 4, "wrong interval for CACA");

  iv = getInterval(esa, "TGT", 3);
  mu_assert(iv.lb == 98 && iv.rb == 100 && iv.lcp == 3, "wrong interval for TGT");

  // only a prefix matched
  iv = getInterval(esa, "GTGTA", 5);
  mu_assert(iv.lb == 89 && iv.rb == 90 && iv.lcp == 4, "wrong interval for GTGTA");

  // nothing matched -> lcp 0, "complete" interval
  iv = getInterval(esa, "X", 1);
  mu_assert(iv.lb == 0 && iv.rb == esa->n - 1 && iv.lcp == 0, "wrong interval for X");

  freeEsa(esa);
  return NULL;
}

// naive: length of lcp of suffixes i and j of given seq (1-indexed)
size_t lcpNaive(char *str, size_t n, size_t i, size_t j) {
  size_t k = 0;
  while (MAX(i, j) + k < n && str[k + i] == str[k + j])
    k++;
  return k;
}

char *test_lcpTreeSmall() {
  size_t n = strlen(seq2);
  Esa *esa = getEsa(seq2, n); // calculate esa, including $
  Interval *tree = getLcpTree(esa);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++) {
      int64_t exp = lcpNaive(seq2, n, esa->sa[i], esa->sa[j]);
      int64_t obs = getLcp(esa, tree, i, j);
      if (exp != obs)
        printf("tried %zu and %zu", i, j);
      mu_assert_eq(exp, obs, "wrong lcp value");
    }

  freeLcpTree(tree);
  freeEsa(esa);
  return NULL;
}

char *test_lcpTree() {
  size_t n = 1000000;
  char *s = randSeq(n);
  Esa *esa = getEsa(s, n + 1);
  Interval *tree = getLcpTree(esa);

  for (size_t i = 0; i < 50; i++) {
    size_t a = rand() % (n - 10) + 5;
    size_t b = a + (rand() % 10) - 5;
    size_t exp = lcpNaive(s, n + 1, a, b);
    size_t obs = getLcp(esa, tree, esa->isa[a], esa->isa[b]);
    if (exp != obs)
      fprintf(stderr, "tried %zu and %zu...\n", a, b);
    mu_assert_eq(exp, obs, "lcp not correct");
  }

  freeLcpTree(tree);
  freeEsa(esa);
  free(s);
  return NULL;
}

char *all_tests() {
  srand(time(NULL));
  mu_suite_start();
  mu_run_test(test_intervals);
  mu_run_test(test_lcpTreeSmall);
  mu_run_test(test_lcpTree);
  return NULL;
}
RUN_TESTS(all_tests)
