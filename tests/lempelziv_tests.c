#include "minunit.h"
#include <string.h>

#include "esa.h"
#include "lempelziv.h"

// hotspot paper example
static char *seq = "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"
                   "CACATATGCTAACTCTCAGTCTGTGTGTGCA$";
// correct factorization
static char *factors[] = {
    "G",   "C",        "A",
    "C",   "GCACGCAC", "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
    "T",   "AT",       "GC",
    "TA",  "AC",       "T",
    "CTC", "A",        "G",
    "TCT", "GT",       "GTGTG",
    "CA",  "$"};

char *test_LempelZiv() {
  size_t n = strlen(seq);
  Esa *esa = getEsa(seq, n); // calculate esa, including $

  Fact *lzf = computeLZFact(esa);
  mu_assert(lzf->n == 20, "wrong number of LZ factors");

  for (size_t i = 0; i < lzf->n; i++) {
    mu_assert(!strncmp(lzf->str + lzf->fact[i], factors[i], factLen(lzf, i)),
              "wrong factor");
  }

  freeFact(lzf);
  freeEsa(esa);
  return NULL;
}

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_LempelZiv);
  return NULL;
}
RUN_TESTS(all_tests)
