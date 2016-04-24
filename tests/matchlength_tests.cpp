#include "minunit.h"
#include <string.h>

#include "esa.h"
#include "matchlength.h"

// hotspot paper example 1
static char const *seq1 = "CCCCGCTCTCCA$";
// hotspot paper example 2
static char const *seq2 =
    "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC"
    "ACACATATGCTAACTCTCAGTCTGTGTGTGCA$";

static char const *factors1[] = {"CCC", "C", "G", "CTC", "TC", "C", "A", "$"};
static char const *factors2[] = {
    "GCACGCAC", "GCAC", "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
    "TA",       "TGC",  "TA",
    "AC",       "TCT",  "CA",
    "GT",       "CT",   "GTGTG",
    "TGC",      "A",    "$"};

char const *checkML(char const *seq, char const *facts[], size_t num) {
  Esa esa(seq, strlen(seq)); // calculate esa, including $

  Fact mlf = computeMLFact(esa);
  mu_assert(mlf.n == num, "wrong number of ML factors");

  for (size_t i = 0; i < mlf.n; i++) {
    mu_assert(!strncmp(mlf.str + mlf.fact[i], facts[i], factLen(mlf, i)), "wrong factor");
  }

  return NULL;
}

char const *test_MatchLength1() { return checkML(seq1, factors1, 8); }
char const *test_MatchLength2() { return checkML(seq2, factors2, 15); }

char const *all_tests() {
  mu_suite_start();
  mu_run_test(test_MatchLength1);
  mu_run_test(test_MatchLength2);
  return NULL;
}
RUN_TESTS(all_tests)
