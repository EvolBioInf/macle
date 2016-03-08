#include "minunit.h"
#include <string.h>

#include "sequenceData.h"
#include "esa.h"
#include "matchlength.h"

// hotspot paper example 1
static char *seq1 = "CCCCGCTCTCCA$";
// hotspot paper example 2
static char *seq2 = "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC"
                    "ACACATATGCTAACTCTCAGTCTGTGTGTGCA$";

static char *factors1[] = {"CCC", "C", "G", "CTC", "TC", "C", "A", "$"};
static char *factors2[] = {
    "GCACGCAC", "GCAC", "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
    "TA",       "TGC",  "TA",
    "AC",       "TCT",  "CA",
    "GT",       "CT",   "GTGTG",
    "TGC",      "A",    "$"};

char *checkML(char *seq, char *facts[], size_t num) {
  Esa *esa = getEsa(seq, strlen(seq)); // calculate esa, including $

  Fact *mlf = computeMLFact(esa);
  mu_assert(mlf->n == num, "wrong number of ML factors");

  for (size_t i = 0; i < mlf->n; i++) {
    mu_assert(!strncmp(mlf->str + mlf->fact[i], facts[i], factLen(mlf, i)),
              "wrong factor");
  }

  freeFact(mlf);
  freeEsa(esa);
  return NULL;
}

char *test_MatchLength1() { return checkML(seq1, factors1, 8); }
char *test_MatchLength2() { return checkML(seq2, factors2, 15); }

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_MatchLength1);
  mu_run_test(test_MatchLength2);
  return NULL;
}
RUN_TESTS(all_tests)
