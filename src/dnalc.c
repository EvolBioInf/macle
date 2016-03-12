#include "prelude.h"

#include "args.h"
#include "eprintf.h"
#include "gsl_rng.h"

#include "sequenceData.h"
#include "matchlength.h"
#include "lempelziv.h"
#include "list.h"
#include "periodicity.h"
#include "complexity.h"

// ---- for benchmarking ----
#include <sys/resource.h>
struct rusage ruse;
double last_tick_time = 0;
#define CPU_TIME                                                                         \
  (getrusage(RUSAGE_SELF, &ruse),                                                        \
   ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +                                         \
       1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
void tick() { last_tick_time = CPU_TIME; }
void tock(bool b, char *str) {
  if (b)
    fprintf(stderr, "[BENCH] %s: %.2fs\n", str, CPU_TIME - last_tick_time);
}
// ---------------------------

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

void printPlot(uint32_t k, size_t n, Sequence *seq, double **ys) {
  printf("offset ");
  for (int i = 0; i < seq->numSeq; i++) { // two columns for each seq
    printf("\"%s\" ", seq->headers[i] + 1);
    printf("\"%s\" ", seq->headers[i] + 1);
  }
  printf("\n");
  for (size_t j = 0; j < n; j++) {
    printf("%zu ", j * k);
    for (int i = 0; i < 2 * seq->numSeq; i++)
      printf("%.4f ", ys[i][j]);
    printf("\n");
  }
}

// read sequence, do stuff
void scanFile(Sequence *seq) {
  bool b = args.b;
  size_t maxlen = maxSeqLen(seq); // this implies the domain of the plot
  size_t minlen = minSeqLen(seq); // this restricts the reasonable window sizes

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
  double **ys = emalloc(2 * seq->numSeq * sizeof(double *));
  for (int i = 0; i < 2 * seq->numSeq; i++)
    ys[i] = ecalloc(entries + 1, sizeof(double));

  // TODO: use gc content of complete file or just current sequence?
  double gc = gcContent(seq);

  for (int i = 0; i < seq->numSeq; i++) {
    char *t = seqStr(seq, i);
    size_t n = seqLen(seq, i);

    tick();
    Esa *esa = getEsa(t, n + 1); // esa for sequence+$
    tock(b, "getEsa");
    if (args.p) {
      printf("%s\n", seq->headers[i]);
      /* printEsa(esa); */
    }

    tick();
    Fact *mlf = computeMLFact(esa);
    tock(b, "computeMLFact");
    if (args.p) {
      printf("ML-Factors (%zu):\n", mlf->n);
      printFact(mlf);
    }

    tick();
    mlComplexity(w, k, ys[2 * i], mlf, gc);
    tock(b, "mlComplexity");
    freeFact(mlf);

    tick();
    Fact *lzf = computeLZFact(esa, false);
    tock(b, "computeLZFact");
    tick();
    size_t pnum;
    List **ls = getPeriodicityLists(true, lzf, esa, &pnum);
    tock(b, "getPeriodicityLists");

    tick();
    runComplexity(w, k, ys[2 * i + 1], n + 1, ls);
    tock(b, "runComplexity");

    if (args.p) {
      printf("LZ-Factors (%zu):\n", lzf->n);
      printFact(lzf);
      printf("Periodicities (%zu):\n", pnum);
      for (size_t j = 0; j < n + 1; j++)
        for (eachListItem(curr, ls[j]))
          printPeriodicity((Periodicity *)curr->value);
    }

    freeFact(lzf);
    freePeriodicityLists(ls, n + 1);

    freeEsa(esa);
  }

  if (!args.p) {
    if (args.g) { // print to be directly piped into gnuplot
      gnuplotCode(w, k, 2 * seq->numSeq);
      for (int i = 0; i < 2 * seq->numSeq; i++) {
        printPlot(k, entries, seq, ys);
        printf("e\n");
      }
    } else { // just print resulting data
      printPlot(k, entries, seq, ys);
    }
  }

  for (int i = 0; i < 2 * seq->numSeq; i++)
    free(ys[i]);
  free(ys);
}

int main(int argc, char *argv[]) {
  setprogname2(PROGNAME);
  parseArgs(argc, argv);
  gsl_rng *rng = ini_gsl_rng(args.s); // init seed, if provided

  // process files (or stdin, if none given)
  if (!args.num_files) {
    Sequence *seq = readFastaFromFile(NULL);
    scanFile(seq);
    freeSequence(seq);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      Sequence *seq = readFastaFromFile(args.files[i]);
      scanFile(seq);
      freeSequence(seq);
    }
  }

  free_gsl_rng(rng, args.s); // save seed, if not using provided
}
