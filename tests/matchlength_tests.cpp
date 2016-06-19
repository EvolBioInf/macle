#include "minunit.h"
#include <string>
using namespace std;

#include "esa.h"
#include "matchlength.h"
#include "util.h"

// hotspot paper example 1
static string seq1 = "CCCCGCTCTCCA";
// hotspot paper example 2
static string seq2 =
    "GCACGCACGCACACACACACACACATATGCTAACTCTCAGTCTGTGTGTGCA";

static string factors1[] = {"CCC", "CG", "CTC", "TC", "C", "A", "$"};
static string factors2[] = {
    "GCACGCAC", "GCACACACA", "CACACACA",
    "TATG",     "CT",        "A",
    "ACT",      "CTC",       "AGT",
    "CTG",      "TGTGTGC",   "A",
    "$"};

void checkML(string seq, string facts[], size_t num) {
  string s = seq + "$" + revComp(seq) + "$";
  Esa esa(s.c_str(), s.size()); // calculate esa, both strands
  Fact mlf;
  computeMLFact(mlf, esa);
  cout << s << endl;
  mlf.print();
  mu_assert(mlf.fact.size() == num, "wrong number of ML factors");

  for (size_t i = 0; i < mlf.fact.size(); i++)
    mu_assert_eq(facts[i], string(mlf.str + mlf.fact[i], mlf.factLen(i)), "wrong factor");
}

void test_MatchLength1() { return checkML(seq1, factors1, 7); }
void test_MatchLength2() { return checkML(seq2, factors2, 13); }

void all_tests() {
  mu_run_test(test_MatchLength1);
  mu_run_test(test_MatchLength2);
}
RUN_TESTS(all_tests)
