#include "minunit.h"
#include <string>
using namespace std;

#include "index.h"

void assert_dataEqual(vector<ComplexityData> const &d1, vector<ComplexityData> const &d2, bool onlyInfo) {
  mu_assert_eq(d1.size(), d2.size(), "Number of sequences not equal");
  for (size_t i=0; i<d1.size(); i++) {
    ComplexityData const &c1 = d1[i];
    ComplexityData const &c2 = d2[i];
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
      continue;

    mu_assert_eq(c1.bad.size(), c2.bad.size(), "Num. of bad intervals not equal");
    for (size_t j=0; j<c1.bad.size(); j++) {
      mu_assert_eq(c1.bad[j].first, c2.bad[j].first, "bad interval not equal");
      mu_assert_eq(c1.bad[j].second, c2.bad[j].second, "bad interval not equal");
    }

    mu_assert_eq(c1.mlf.size(), c2.mlf.size(), "Num. of MLF not equal");
    for (size_t j=0; j<c1.mlf.size(); j++)
      mu_assert_eq(c1.mlf[j], c2.mlf[j], "MLF not equal");

    mu_assert_eq(c1.pl.size(), c2.pl.size(), "Num. of PL not equal");
    for (size_t j=0; j<c1.pl.size(); j++)
      for (size_t k=0; k<c1.pl[j].size(); k++) {
        mu_assert_eq(c1.pl[j][k].b, c2.pl[j][k].b, "Per (b) not equal");
        mu_assert_eq(c1.pl[j][k].e, c2.pl[j][k].e, "Per (e) not equal");
        mu_assert_eq(c1.pl[j][k].l, c2.pl[j][k].l, "Per (l) not equal");
      }
  }
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

  //create 2 index files - joined and unjoined
  vector<ComplexityData> dat;
  extractData(dat,ff,false);
  mu_assert_eq((size_t)2, dat.size(), "wrong number of sequences");
  mu_assert_eq((size_t)0, dat[0].regions.size(), "wrong number of regions");
  mu_assert_eq((size_t)0, dat[1].regions.size(), "wrong number of regions");
  mu_assert_eq(seq1.name, dat[0].name, "wrong sequence name");
  mu_assert_eq(seq2.name, dat[1].name, "wrong sequence name");
  saveData(dat, iname1);

  vector<ComplexityData> datJ;
  extractData(datJ,ff,true);
  mu_assert_eq((size_t)1, datJ.size(), "wrong number of sequences");
  mu_assert_eq((size_t)2, datJ[0].regions.size(), "wrong number of regions");
  mu_assert_eq(ff.filename, datJ[0].name, "wrong sequence name");
  saveData(datJ, iname2);

  //try loading and check
  cerr << "load unjoined complete..." << endl;
  vector<ComplexityData> dat2;
  loadData(dat2, iname1, false);
  assert_dataEqual(dat, dat2, false);

  cerr << "load joined complete..." << endl;
  vector<ComplexityData> datJ2;
  loadData(datJ2, iname2, false);
  assert_dataEqual(datJ, datJ2, false);

  //try loading only info and check
  cerr << "load unjoined info..." << endl;
  dat2.clear();
  loadData(dat2, iname1, true);
  assert_dataEqual(dat, dat2, true);

  cerr << "load joined info..." << endl;
  datJ2.clear();
  loadData(datJ2, iname2, true);
  assert_dataEqual(datJ, datJ2, true);

  remove(iname1);
  remove(iname2);
}

void all_tests() {
  mu_run_test(test_saveLoadData);
}
RUN_TESTS(all_tests)
