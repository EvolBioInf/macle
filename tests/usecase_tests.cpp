#include "minunit.h"
#include <string>
using namespace std;

#include "index.h"
#include "complexity.h"

ResultMat ys;
FastaFile ff;
ComplexityData datJ;

// dnalc seq.fa, dnalc -j seq.fa
void test_global_no_settings() {
  size_t w=0, k=0;
  //joined, global, no sequence specified
  ys = calcComplexities(w, k, Task(-1, 0, 0), datJ);
  mu_assert_eq(1UL, ys.size(), "wrong number of columns!");
  mu_assert_eq(1UL, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ.name+" (MC)", ys[0].first, "wrong header!");
}

void test_global_with_settings() {
  size_t w=0, k=0;
  //joined, global, all seqs
  ys = calcComplexities(w, k, Task(-1,0,0), datJ);
  mu_assert_eq(1UL, ys.size(), "wrong number of columns!");
  mu_assert_eq(1UL, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ.name+" (MC)", ys[0].first, "wrong header!");

  //joined, global, second sequence
  ys = calcComplexities(w, k, Task(1,0,0), datJ);
  mu_assert_eq(1UL, ys.size(), "wrong number of columns!");
  mu_assert_eq(1UL, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ.labels[1]+" (MC)", ys[0].first, "wrong header!");
}

void all_tests() {
  //construct a sequence file
  FastaSeq seq1("seq1","comment","NNNNNATATATGCGCGCATGCATGCNNNNN");
  FastaSeq seq2("seq2","comment","NNNNNNNNNNNATCGACATGCTANNNNGTGAGTCTANNNN");
  ff.filename = "seq.fa";
  ff.seqs.push_back(seq1);
  ff.seqs.push_back(seq2);
  extractData(datJ,ff);

  mu_run_test(test_global_no_settings);
  mu_run_test(test_global_with_settings);
}
RUN_TESTS(all_tests)
