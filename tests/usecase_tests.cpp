#include "minunit.h"
#include <string>
using namespace std;

#include "index.h"
#include "complexity.h"

ResultMat ys;
FastaFile ff;
vector<ComplexityData> dat;
vector<ComplexityData> datJ;

// dnalc seq.fa, dnalc -j seq.fa
void test_global_no_settings() {
  size_t w=0, k=0;
  //unjoined, global, both compl., no sequence specified
  ys = calcComplexities(w, k, 'b', 0, dat);
  mu_assert_eq((size_t)0, w, "w not zero anymore!");
  mu_assert_eq((size_t)0, k, "w not zero anymore!");
  mu_assert_eq((size_t)4, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(dat[0].name+" (MC)", ys[0].first, "wrong header!");
  mu_assert_eq(dat[0].name+" (RC)", ys[1].first, "wrong header!");
  mu_assert_eq(dat[1].name+" (MC)", ys[2].first, "wrong header!");
  mu_assert_eq(dat[1].name+" (RC)", ys[3].first, "wrong header!");

  //joined, global, both compl., no sequence specified
  ys = calcComplexities(w, k, 'b', 0, datJ);
  mu_assert_eq((size_t)2, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ[0].name+" (MC)", ys[0].first, "wrong header!");
  mu_assert_eq(datJ[0].name+" (RC)", ys[1].first, "wrong header!");
}

void test_global_with_settings() {
  size_t w=0, k=0;
  ys = calcComplexities(w, k, 'b', 0, dat);
  //unjoined, global, both compl., second sequence
  auto ys2 = calcComplexities(w, k, 'b', 2, dat);
  mu_assert_eq((size_t)2, ys2.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys2[0].second.size(), "wrong number of entries!");
  mu_assert_eq(ys2[0].first, dat[1].name+" (MC)", "wrong header!");
  mu_assert_eq(ys2[1].first, dat[1].name+" (RC)", "wrong header!");
  mu_assert_eq(ys[2].first, ys2[0].first, "wrong column content!");
  mu_assert_eq(ys[3].first, ys2[1].first, "wrong column content!");
  mu_assert_eq(ys[2].second[0], ys2[0].second[0], "wrong column content!");
  mu_assert_eq(ys[3].second[0], ys2[1].second[0], "wrong column content!");

  //joined, global, both compl., second sequence
  ys = calcComplexities(w, k, 'b', 2, datJ);
  mu_assert_eq((size_t)2, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ[0].labels[1]+" (MC)", ys[0].first, "wrong header!");
  mu_assert_eq(datJ[0].labels[1]+" (RC)", ys[1].first, "wrong header!");

  //unjoined, global, all seqs, just RC
  ys = calcComplexities(w, k, 'r', 0, dat);
  mu_assert_eq((size_t)2, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(dat[0].name+" (RC)", ys[0].first, "wrong header!");
  mu_assert_eq(dat[1].name+" (RC)", ys[1].first, "wrong header!");

  //joined, global, all seqs, just MC
  ys = calcComplexities(w, k, 'm', 0, datJ);
  mu_assert_eq((size_t)1, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ[0].name+" (MC)", ys[0].first, "wrong header!");

  //joined, global, second sequence, just RC
  ys = calcComplexities(w, k, 'r', 2, datJ);
  mu_assert_eq((size_t)1, ys.size(), "wrong number of columns!");
  mu_assert_eq((size_t)1, ys[0].second.size(), "wrong number of entries!");
  mu_assert_eq(datJ[0].labels[1]+" (RC)", ys[0].first, "wrong header!");
}

void all_tests() {
  //construct a sequence file
  FastaSeq seq1("seq1","comment","NNNNNATATATGCGCGCATGCATGCNNNNN");
  FastaSeq seq2("seq2","comment","NNNNNNNNNNNATCGACATGCTANNNNGTGAGTCTANNNN");
  ff.filename = "seq.fa";
  ff.seqs.push_back(seq1);
  ff.seqs.push_back(seq2);

  //get unjoined and joined data
  extractData(dat,ff,false);
  extractData(datJ,ff,true);

  mu_run_test(test_global_no_settings);
  mu_run_test(test_global_with_settings);
}
RUN_TESTS(all_tests)
