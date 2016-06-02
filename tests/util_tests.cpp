#include "minunit.h"
#include <vector>
using namespace std;

#include "util.h"

void test_randSeq() {
  size_t n = 80;
  string alphabet = "XYZ";
  string seq = randSeq(n, alphabet);
  mu_assert_eq(n, seq.size(), "incorrect size!");
  cout << seq << endl;
  for (char c : seq)
    mu_assert(alphabet.find(c) != string::npos, "invalid character: " << c);
  // test with default alphabet, different gc contents
  alphabet="ACGT";
  seq = randSeq(n, 0.4);
  cout << seq << endl;
  for (char c : seq)
    mu_assert(alphabet.find(c) != string::npos, "invalid character: " << c);
  alphabet="AT";
  seq = randSeq(n, 0.0);
  cout << seq << endl;
  for (char c : seq)
    mu_assert(alphabet.find(c) != string::npos, "invalid character: " << c);
  alphabet="GC";
  seq = randSeq(n, 1.0);
  cout << seq << endl;
  for (char c : seq)
    mu_assert(alphabet.find(c) != string::npos, "invalid character: " << c);
}

void test_randRun() {
  size_t n = 1000;
  for (auto l : vector<size_t>{1, 2, 4, 8, 16}) {
    string runseq = randRun(n, l);
    string per = runseq.substr(0, l);
    size_t cnt = 0;
    size_t pos = 0;
    while (runseq.find(per, pos) != string::npos) {
      cnt++;
      pos += l;
    }
    mu_assert_eq(n / l, cnt, "wrong number of repeats");
  }
}

// test reverse complement. assumes normalized sequence, so just ACGT toggle
// and whole string is reversed.
void test_revComp() {
  string test = "ACCCGGGGGNNNNNNNTTTTTTTatcg$";
  string result = revComp(test);
  mu_assert(result == "$gctaAAAAAAANNNNNNNCCCCCGGGT", "wrong reverse complement!");
}

void all_tests() {
  mu_run_test(test_randSeq);
  mu_run_test(test_randRun);
  mu_run_test(test_revComp);
}
RUN_TESTS(all_tests)
