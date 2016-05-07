#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

#include "args.h"
#include "bench.h"

#include "fastafile.h"
#include "matchlength.h"
#include "lempelziv.h"
#include "periodicity.h"
#include "complexity.h"
#include "util.h"

using namespace std;

// All information from a sequence required to calculate complexity plots
struct ComplexityData {
  void print() const;
  std::string name;                       // name of sequence
  size_t len;                             // length of sequence
  double gc;                              // gc content of sequence
  std::vector<size_t> mlf;                // match factors
  std::vector<std::list<Periodicity>> pl; // periodicities
};

void ComplexityData::print() const {
  cout << name << " " << len << " " << gc << endl;
  cout << mlf.size() << endl;
  for (auto i : mlf)
    cout << i << endl;
  size_t pernum = 0;
  for (auto &l : pl)
    pernum += l.size();
  cout << pernum << endl;
  for (auto &l : pl)
    for (auto p : l)
      cout << p.b << " " << p.e << " " << p.l << endl;
}

bool loadData(vector<ComplexityData> &cplx, char const *file) {
  istream *finP = &cin;
  ifstream fs;
  if (file) {
    fs = ifstream(file);
    if (!fs.is_open()) {
      cerr << "ERROR: Could not open file: " << file << endl;
      return false;
    }
    finP = &fs;
  }
  istream &fin = *finP;

  int n;
  fin >> n;
  for (int i = 0; i < n; i++) {
    ComplexityData c;
    fin >> c.name >> c.len >> c.gc;
    int fnum;
    fin >> fnum;
    c.mlf.resize(fnum);
    for (int j = 0; j < fnum; j++)
      fin >> c.mlf[j];
    int pnum;
    fin >> pnum;
    c.pl.resize(c.len);
    for (int j = 0; j < pnum; j++) {
      size_t b, e, l;
      fin >> b >> e >> l;
      c.pl[b].push_back(Periodicity(b, e, l));
    }
    cplx.push_back(c);
  }
  return true;
}

// given sequences from a fasta file, calculate match factors and runs
void extractData(vector<ComplexityData> &cplx, FastaFile &ff) {
  for (auto &seq : ff.seqs) {
    ComplexityData c;
    c.name = seq.name;
    c.gc = gcContent(seq.seq);
    c.len = seq.seq.size();

    string &s = seq.seq;
    s = s + "$" + revComp(s) + "$";
    tick();
    Esa esa(s.c_str(), s.size()); // esa for seq+$+revseq+$
    tock("getEsa (2n)");
    tick();
    Fact mlf;
    computeMLFact(mlf, esa);
    tock("computeMLFact");
    c.mlf = mlf.fact;

    if (args.p) {
      cout << "ML-Factors (" << mlf.fact.size() << "):" << endl;
      mlf.print();
    }

    tick();
    s.resize(s.size() / 2); // drop complementary seq.
    s.shrink_to_fit();
    esa.str = s.c_str();
    esa.n = s.size();
    reduceEsa(esa);
    tock("reduceEsa");

    if (args.p) {
      cout << seq.name << seq.comment << endl;
      // esa.print();
    }

    tick();
    Fact lzf;
    computeLZFact(lzf, esa, false);
    tock("computeLZFact");
    tick();
    size_t pnum;
    c.pl = getPeriodicityLists(true, lzf, esa, pnum);
    tock("getPeriodicityLists");

    if (args.p) {
      cout << "LZ-Factors (" << lzf.fact.size() << "):" << endl;
      lzf.print();
      cout << "Periodicities (" << pnum << "):" << endl;
      for (auto &l : c.pl)
        for (auto p : l)
          printPeriodicity(p);
    }

    cplx.push_back(c);
  }
}

void gnuplotCode(uint32_t w, uint32_t k, int n) {
  printf( //"set terminal png; "
      "set key autotitle columnheader; "
      "set ylabel \"window complexity\"; "
      "set xlabel \"window offset (w=%u, k=%u)\"; ",
      w, k);
  for (int i = 0; i < n; i++)
    printf("%s using 1:%d with lines", i ? ", ''" : "plot \"-\"", i + 2);
  printf(";\n");
}

// print data: X Y1 ... Yn
void printPlot(uint32_t w, uint32_t k, size_t n, vector<vector<double>> &ys) {
  for (size_t j = 0; j < n; j++) {
    printf("%zu\t", j * k + w / 2); // center of window
    for (size_t i = 0; i < ys.size(); i++)
      printf("%.4f\t", ys[i][j]);
    printf("\n");
  }
}

// read data, do stuff
void processFile(char const *file) {
  vector<ComplexityData> dat;
  if (args.l) {
    tick();
    if (!loadData(dat, file))
      return;
    tock("loadData");
  } else { // not loading from pre-computed data -> fasta file
    tick();
    FastaFile ff(file);
    tock("readFastaFromFile");
    if (ff.failed) {
      cerr << "Skipping invalid FASTA file..." << endl;
      return;
    }
    extractData(dat, ff);

    if (args.s && !args.p) { // just dump intermediate results and quit
      cout << dat.size() << endl;
      for (auto d : dat)
        d.print();
      return;
    }
  }

  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  for (size_t i = 0; i < dat.size(); i++) {
    size_t len = dat[i].len;
    if (len > maxlen)
      maxlen = len;
    if (len < minlen)
      minlen = len;
  }

  // adapt window size and interval
  size_t w = args.w;
  size_t k = args.k;
  if (w == 0)
    w = minlen;       // default window = whole (smallest) seq.
  w = min(w, minlen); // biggest window = whole (smallest) seq.
  if (k == 0)
    k = w / 10;  // default interval = 1/10 of window
  k = min(k, w); // biggest interval = window size

  // array for results for all sequences in file
  size_t entries = (maxlen - w) / k + 1;
  vector<vector<double>> ys(2 * dat.size(), vector<double>(entries + 1));

  int n = 0;
  for (auto &seq : dat) {
    tick();
    mlComplexity(w, k, ys[2 * n], seq.len, seq.mlf, seq.gc);
    tock("mlComplexity");

    tick();
    runComplexity(w, k, ys[2 * n + 1], seq.len, seq.pl);
    tock("runComplexity");

    n++;
  }

  if (!args.p) {
    if (args.g) { // print to be directly piped into some plot tool
      if (args.gf == 2) {
        for (size_t i = 0; i < 2 * dat.size(); i++) {
          for (size_t j = 0; j < entries; j++)
            printf("%zu\t%.4f\n", j * k + w / 2, ys[i][j]);
          printf("\n");
        }
      } else if (args.gf == 1) {
        gnuplotCode(w, k, 2 * dat.size()); // gnuplot control code
        // need to output everything n times for n plots
        for (size_t i = 0; i < 2 * dat.size(); i++) {
          // print column header (for plot labels)
          printf("offset ");
          for (size_t j = 0; j < dat.size(); j++) { // two columns for each seq
            printf("\"%s %s\"\t", dat[j].name.c_str(), "(MC)");
            printf("\"%s %s\"\t", dat[j].name.c_str(), "(PC)");
          }
          printf("\n");
          printPlot(w, k, entries, ys); // print plot itself
          printf("e\n");                // separator between repeats of data
        }
      }
    } else { // just print resulting data
      printPlot(w, k, entries, ys);
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  // process files (or stdin, if none given)
  tick();
  if (!args.num_files) {
    processFile(nullptr);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      processFile(args.files[i]);
    }
  }
  tock("total time");
}
