#include "minunit.h"
#include <string>
using namespace std;

#include "index.h"

void assert_dataEqual(ComplexityData const &c1, ComplexityData const &c2, bool onlyInfo) {
  mu_assert_eq(c1.name, c2.name, "Names not equal!");
  mu_assert_eq(c1.len, c2.len, "Length not equal!");
  mu_assert_eq(c1.gc, c2.gc, "GC not equal!");

  mu_assert_eq(c1.labels.size(), c1.regions.size(), "Num. Labels and Regions different");
  mu_assert_eq(c2.labels.size(), c2.regions.size(), "Num. Labels and Regions different");

  mu_assert_eq(c1.labels.size(), c2.labels.size(), "Num. of labels not equal");
  mu_assert_eq(c1.regions.size(), c2.regions.size(), "Num. of labels not equal");
  for (size_t j=0; j<c1.labels.size(); j++) {
    mu_assert_eq(c1.labels[j], c2.labels[j], "Label not equal");
    mu_assert_eq(c1.regions[j].first, c2.regions[j].first, "Region not equal");
    mu_assert_eq(c1.regions[j].second, c2.regions[j].second, "Region not equal");
  }
  mu_assert_eq(c1.numbad, c2.numbad, "Numbad not equal");
  if (onlyInfo)
    return;

  mu_assert_eq(c1.bad.size(), c2.bad.size(), "Num. of bad intervals not equal");
  for (size_t j=0; j<c1.bad.size(); j++) {
    mu_assert_eq(c1.bad[j].first, c2.bad[j].first, "bad interval not equal");
    mu_assert_eq(c1.bad[j].second, c2.bad[j].second, "bad interval not equal");
  }

  mu_assert_eq(c1.mlf.size(), c2.mlf.size(), "Num. of MLF not equal");
  for (size_t j=0; j<c1.mlf.size(); j++)
    mu_assert_eq(c1.mlf[j], c2.mlf[j], "MLF not equal");
}

void test_saveLoadData() {
  char const* iname1 = "_tmp_seq.fa.bin";
  char const* iname2 = "_tmp_seq.fa.joined.bin";

  FastaSeq seq1("seq1","comment","NNNNNATATATGCGCGCATGCATGCNNNNN");
  FastaSeq seq2("seq2","comment","NNNNNNNNNNNATCGACATGCTANNNNGTGAGTCTANNNN");
  FastaFile ff;
  ff.filename = "seq.fa";
  ff.seqs.push_back(seq1);
  ff.seqs.push_back(seq2);

  //create index file
  ComplexityData datJ;
  extractData(datJ,ff);
  mu_assert_eq((size_t)2, datJ.regions.size(), "wrong number of regions");
  mu_assert_eq(ff.filename, datJ.name, "wrong sequence name");
  saveData(datJ, iname2);

  //try loading and check
  cerr << "load joined complete..." << endl;
  ComplexityData datJ2;
  loadData(datJ2, iname2, false);
  assert_dataEqual(datJ, datJ2, false);

  //try loading only info and check
  cerr << "load joined info..." << endl;
  datJ2 = ComplexityData();
  loadData(datJ2, iname2, true);
  assert_dataEqual(datJ, datJ2, true);

  remove(iname1);
  remove(iname2);
}

void all_tests() {
  mu_run_test(test_saveLoadData);
}
RUN_TESTS(all_tests)
