#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"

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

char *all_tests() {
  mu_suite_start();

  mu_run_test(test_getEsa);

  return NULL;
}

RUN_TESTS(all_tests)
