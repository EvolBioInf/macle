#include "minunit.h"
#include <cstring>
using namespace std;

#include "esa.h"
#include "matchlength.h"

// hotspot paper example 1
static string seq1 = "CCCCGCTCTCCA$";
// hotspot paper example 2
static string seq2 =
    "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC"
    "ACACATATGCTAACTCTCAGTCTGTGTGTGCA$";

static string factors1[] = {"CCC", "C", "G", "CTC", "TC", "C", "A", "$"};
static string factors2[] = {
    "GCACGCAC", "GCAC", "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
    "TA",       "TGC",  "TA",
    "AC",       "TCT",  "CA",
    "GT",       "CT",   "GTGTG",
    "TGC",      "A",    "$"};

void checkML(string seq, string facts[], size_t num) {
  Esa esa(seq.c_str(), seq.size()); // calculate esa, including $

  Fact mlf;
  computeMLFact(mlf, esa);
  mu_assert(mlf.fact.size() == num, "wrong number of ML factors");

  for (size_t i = 0; i < mlf.fact.size(); i++) {
    mu_assert(!strncmp(mlf.str + mlf.fact[i], facts[i].c_str(), factLen(mlf, i)),
              "wrong factor");
  }
}

void test_MatchLength1() { return checkML(seq1, factors1, 8); }
void test_MatchLength2() { return checkML(seq2, factors2, 15); }

void all_tests() {
  mu_run_test(test_MatchLength1);
  mu_run_test(test_MatchLength2);
}
RUN_TESTS(all_tests)
