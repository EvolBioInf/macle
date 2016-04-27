#include "args.h"
#include "bench.h"

#include "fastafile.h"
#include "matchlength.h"
#include "lempelziv.h"
#include "periodicity.h"
#include "complexity.h"

using namespace std;

// calculate the GC content
double gcContent(pfasta_seq const &ps) {
  size_t gc = 0;
  char const *ptr = ps.seq;
  for (; *ptr; ptr++)
    if (*ptr == 'g' || *ptr == 'G' || *ptr == 'c' || *ptr == 'C')
      gc++;
  return (double)gc / (ptr - ps.seq);
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

void printPlot(uint32_t k, size_t n, FastaFile &ff, vector<vector<double>> &ys) {
  printf("offset ");
  for (size_t i = 0; i < ff.seqs.size(); i++) { // two columns for each seq
    printf("\"%s %s\" ", ff.seqs[i].name, "(ML)");
    printf("\"%s %s\" ", ff.seqs[i].name, "(Per)");
  }
  printf("\n");
  for (size_t j = 0; j < n; j++) {
    printf("%zu ", j * k);
    for (size_t i = 0; i < 2 * ff.seqs.size(); i++)
      printf("%.4f ", ys[i][j]);
    printf("\n");
  }
}

// read sequence, do stuff
void scanFile(FastaFile &ff) {
  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  for (size_t i = 0; i < ff.seqs.size(); i++) {
    size_t len = ff.seqs[i].len;
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
  w = MIN(w, minlen); // biggest window = whole (smallest) seq.
  if (k == 0)
    k = w;       // default interval = whole (smallest) seq.
  k = MIN(k, w); // biggest interval = window size

  // array for results for all sequences in file
  size_t entries = (maxlen - w) / k + 1;
  vector<vector<double>> ys(2 * ff.seqs.size(), vector<double>(entries + 1));

  for (size_t i = 0; i < ff.seqs.size(); i++) {
    double gc = gcContent(ff.seqs[i]);
    char *t = ff.seqs[i].seq;
    size_t n = ff.seqs[i].len;

    tick();
    Esa esa(t, n + 1); // esa for sequence+$
    tock("getEsa");
    if (args.p) {
      printf("%s %s\n", ff.seqs[i].name, ff.seqs[i].comment);
      /* printEsa(esa); */
    }

    tick();
    Fact mlf = computeMLFact(esa);
    tock("computeMLFact");
    if (args.p) {
      printf("ML-Factors (%zu):\n", mlf.n);
      mlf.print();
    }

    tick();
    mlComplexity(w, k, ys[2 * i], mlf, gc);
    tock("mlComplexity");

    tick();
    Fact lzf = computeLZFact(esa, false);
    tock("computeLZFact");
    tick();
    size_t pnum;
    vector<list<Periodicity>> ls = getPeriodicityLists(true, lzf, esa, pnum);
    tock("getPeriodicityLists");

    tick();
    runComplexity(w, k, ys[2 * i + 1], esa.n, ls);
    tock("runComplexity");

    if (args.p) {
      printf("LZ-Factors (%zu):\n", lzf.n);
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
            printf("%zu %.4f\n", j * k, ys[i][j]);
          printf("\n");
        }
      } else if (args.gf == 1) {
        gnuplotCode(w, k, 2 * ff.seqs.size());
        for (size_t i = 0; i < 2 * ff.seqs.size(); i++) {
          printPlot(k, entries, ff, ys);
          printf("e\n");
        }
      }
    } else { // just print resulting data
      printPlot(k, entries, ff, ys);
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  // process files (or stdin, if none given)
  if (!args.num_files) {
    FastaFile ff(nullptr);
    scanFile(ff);
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
      scanFile(ff);
      tock("total for file");
    }
  }
}
