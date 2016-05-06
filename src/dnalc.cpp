#include <cstdio>
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

void printPlot(uint32_t w, uint32_t k, size_t n, FastaFile &ff,
               vector<vector<double>> &ys) {
  for (size_t j = 0; j < n; j++) {
    printf("%zu\t", j * k + w / 2); // center of window
    for (size_t i = 0; i < 2 * ff.seqs.size(); i++)
      printf("%.4f\t", ys[i][j]);
    printf("\n");
  }
}

// read sequence, do stuff
void processFile(FastaFile &ff) {
  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  for (size_t i = 0; i < ff.seqs.size(); i++) {
    size_t len = ff.seqs[i].seq.size();
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
  vector<vector<double>> ys(2 * ff.seqs.size(), vector<double>(entries + 1));

  for (size_t i = 0; i < ff.seqs.size(); i++) {
    double gc = gcContent(ff.seqs[i].seq);
    string &s = ff.seqs[i].seq;
    s = s + "$" + revComp(s) + "$";

    tick();
    Esa esa(s.c_str(), s.size()); // esa for seq+$+revseq+$
    tock("getEsa (2n)");
    tick();
    Fact mlf;
    computeMLFact(mlf, esa);
    tock("computeMLFact");

    if (args.p) {
      printf("ML-Factors (%zu):\n", mlf.fact.size());
      mlf.print();
    }

    tick();
    mlComplexity(w, k, ys[2 * i], mlf, gc);
    tock("mlComplexity");

    tick();
    s.resize(s.size()/2); //drop complementary seq.
    s.shrink_to_fit();
    esa.str = s.c_str();
    esa.n = s.size();
    reduceEsa(esa);
    tock("reduceEsa");

    if (args.p) {
      printf("%s %s\n", ff.seqs[i].name.c_str(), ff.seqs[i].comment.c_str());
      esa.print();
    }

    tick();
    Fact lzf;
    computeLZFact(lzf, esa, false);
    tock("computeLZFact");
    tick();
    size_t pnum;
    vector<list<Periodicity>> ls = getPeriodicityLists(true, lzf, esa, pnum);
    tock("getPeriodicityLists");

    tick();
    runComplexity(w, k, ys[2 * i + 1], esa.n, ls);
    tock("runComplexity");

    if (args.p) {
      printf("LZ-Factors (%zu):\n", lzf.fact.size());
      lzf.print();
      printf("Periodicities (%zu):\n", pnum);
      for (size_t j = 0; j < ls.size(); j++)
        for (auto it = ls[j].begin(); it != ls[j].end(); it++)
          printPeriodicity(*it);
    }
  }

  if (!args.p) {
    if (args.g) { // print to be directly piped into gnuplot
      if (args.gf == 0) {
        for (size_t i = 0; i < 2 * ff.seqs.size(); i++) {
          for (size_t j = 0; j < entries; j++)
            printf("%zu\t%.4f\n", j * k + w / 2, ys[i][j]);
          printf("\n");
        }
      } else if (args.gf == 1) {
        gnuplotCode(w, k, 2 * ff.seqs.size()); // gnuplot control code
        // need to output everything n times for n plots
        for (size_t i = 0; i < 2 * ff.seqs.size(); i++) {
          // print column header (for plot labels)
          printf("offset ");
          for (size_t j = 0; j < ff.seqs.size(); j++) { // two columns for each seq
            printf("\"%s %s\"\t", ff.seqs[j].name.c_str(), "(MC)");
            printf("\"%s %s\"\t", ff.seqs[j].name.c_str(), "(PC)");
          }
          printf("\n");
          printPlot(w, k, entries, ff, ys); // print plot itself
          printf("e\n");                    // separator between repeats of data
        }
      }
    } else { // just print resulting data
      printPlot(w, k, entries, ff, ys);
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  // process files (or stdin, if none given)
  if (!args.num_files) {
    FastaFile ff(nullptr);
    processFile(ff);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      tick();
      tick();
      FastaFile ff(args.files[i]);
      tock("readFastaFromFile");
      if (ff.failed) {
        fprintf(stderr, "Skipping invalid FASTA file...\n");
        continue;
      }
      processFile(ff);
      tock("total for file");
    }
  }
}
