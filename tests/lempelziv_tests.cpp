#include "minunit.h"
#include "util.h"
#include <cinttypes>
#include <cstring>
#include <ctime>

#include <string>
#include <vector>
using namespace std;

#include "esa.h"
#include "lempelziv.h"

// hotspot paper example
string seq = "GCACGCACGCACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA"
             "CACATATGCTAACTCTCAGTCTGTGTGTGCA$";

// correct factorization
string factors[] = {
    "G",   "C",        "A",
    "C",   "GCACGCAC", "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
    "T",   "AT",       "GC",
    "TA",  "AC",       "T",
    "CTC", "A",        "G",
    "TCT", "GT",       "GTGTG",
    "CA",  "$"};

// prevOcc example from Crochemore/Ilie paper
string str = "abbaabbbaaabab";
/* static int64_t poResult[] = {3,3,0,10,0,-1,0,7,2,1,2,1,-1,1}; //prevOcc[sa[i]] */

void test_LempelZiv() {
  Esa esa(seq.c_str(), seq.size()); // calculate esa, including $
  auto oldlcp = esa.lcp;

  Fact lzf;
  computeLZFact(lzf, esa, false);
  for (size_t i = 0; i < esa.n + 1; i++) // this was actually a bug!
    mu_assert_eq(oldlcp[i], esa.lcp[i], "old lcp values were changed!");

  for (size_t i = 0; i < lzf.fact.size(); i++) {
    mu_assert(!strncmp(lzf.str + lzf.fact[i], factors[i].c_str(), factLen(lzf, i)),
              "wrong factor");
  }

  mu_assert_eq((size_t)20, lzf.fact.size(), "wrong number of LZ factors");
}

void test_prevOcc() {
  Esa esa(str.c_str(), str.size());
  Fact lzf;
  computeLZFact(lzf, esa, false);
  Fact lzfRef;
  computeLZFact(lzfRef, esa, true);

  // should be same as poResult[sa[i]]
  for (size_t i = 0; i < str.size(); i++)
    mu_assert(lzf.prevOcc[i] == lzfRef.prevOcc[i], "wrong prevOcc array");
}

void test_randomSequence() {
  size_t n = 1000000;
  string ss = randSeq(n);
  ss += "$";
  char const *s = ss.c_str();
  fprintnf(stdout, s, 80);
  printf("\n");
  Esa esa(s, ss.size());
  Fact lzf;
  computeLZFact(lzf, esa, false);
  Fact lzfRef;
  computeLZFact(lzfRef, esa, true); // using different algorithm

  for (size_t i = 0; i < esa.n; i++)
    mu_assert(lzf.prevOcc[i] == lzf.prevOcc[i], "incorrect prevOcc array");
  for (size_t i = 0; i < esa.n; i++)
    mu_assert(lzf.lpf[i] == lzf.lpf[i], "incorrect lpf array");
  mu_assert(lzf.fact.size() == lzfRef.fact.size(), "different number of factors");
  for (size_t i = 0; i < lzf.fact.size(); i++)
    mu_assert(lzf.fact[i] == lzf.fact[i], "incorrect factor positions");

  for (size_t i = 0; i < lzf.fact.size(); i++) {
    int64_t po = lzf.prevOcc[lzf.fact[i]];
    bool tmp = po <= (int64_t)lzf.fact[i];
    if (!tmp)
      printf("prev: %ld curr: %ld\n", po, lzf.fact[i]);
    mu_assert(tmp, "prev can not be after");
    if (po >= 0)
      mu_assert(!strncmp(esa.str + lzf.fact[i], esa.str + po, factLen(lzf, i)),
                "prev is not a match");
  }
}

void all_tests() {
  srand(time(NULL));
  mu_run_test(test_LempelZiv);
  mu_run_test(test_prevOcc);
  for (size_t i = 0; i < 5; i++)
    mu_run_test(test_randomSequence);
}
RUN_TESTS(all_tests)
