#include "prelude.h"

#include <unistd.h>
#include <fcntl.h>

#include "args.h"
#include "eprintf.h"
#include "gsl_rng.h"

#include "sequenceData.h"
#include "complexity.h"
#include "periodicity.h"
#include "lempelziv.h"

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
    fprintf(stderr, "[BENCH] %s: %.3fs\n", str, CPU_TIME - last_tick_time);
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

void printPlot(uint32_t w, uint32_t k, size_t n, Sequence *seq, double **ys) {
  printf("i ");
  for (int i = 0; i < seq->numSeq; i++)
    printf("\"%s\" ", seq->headers[i] + 1);
  printf("\n");
  for (size_t j = 0; j * k <= n - w; j++) {
    printf("%zu ", j * k);
    for (int i = 0; i < seq->numSeq; i++)
      printf("%.4f ", ys[i][j]);
    printf("\n");
  }
}

// read sequence, do stuff
void scanFile(int fd) {
  bool b = args.b;
  tick();
  Sequence *seq = readFasta(fd);
  tock(b, "readFasta");

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

  // array for results
  double **ys = malloc(seq->numSeq * sizeof(double *));
  for (int i = 0; i < seq->numSeq; i++)
    ys[i] = calloc(maxlen, sizeof(double));

  // TODO: use gc content of complete file or just current sequence?
  double gc = gcContent(seq);

  for (int i = 0; i < seq->numSeq; i++) {
    char *t = getSeq(seq, i);
    size_t n = seqLen(seq, i);

    tick();
    Esa *esa = getEsa(t, n+1); //esa for sequence+$
    tock(b, "getEsa");
    tick();
    Fact *mlf = mlComplexity(esa, gc);
    tock(b, "mlComplexity");
    tick();
    Fact *lzf = computeLZFact(esa);
    tock(b, "computeLZFact");

    size_t plen;
    /* Periodicity *ps = getPeriodicities(esa->str, esa->n, &plen); */
    Periodicity *ps = getPeriodicities2(lzf, &plen);
    for (size_t j=0; j<plen; j++)
      printPeriodicity(ps+j);

    if (args.p) {
      printf("%s \t(ML=%.4f)\n", seq->headers[i], mlf->cNor);
      printf("ML-Factors:\n");
      printFact(mlf);
      printf("LZ-Factors:\n");
      printFact(lzf);
    }

    // calculate window complexity, TODO: does this make sense?
    for (size_t j = 0; j * k <= n - w; j++) {
      double facs = factorsFromTo(mlf, j * k, MIN(n, j * k + w) - 1);
      ys[i][j] = (facs / w - mlf->cMin) / (mlf->cMax - mlf->cMin);
    }

    freeFact(mlf);
    freeFact(lzf);

    freeEsa(esa);
  }

  if (!args.p) {
    if (args.g) { // print to be directly piped into gnuplot
      gnuplotCode(w, k, seq->numSeq);
      for (int i = 0; i < seq->numSeq; i++) {
        printPlot(w, k, maxlen, seq, ys);
        printf("e\n");
      }
    } else { // just print resulting data
      printPlot(w, k, maxlen, seq, ys);
    }
  }

  for (int i = 0; i < seq->numSeq; i++)
    free(ys[i]);
  free(ys);
  freeSequence(seq);
}

int main(int argc, char *argv[]) {
  setprogname2(PROGNAME);
  parseArgs(argc, argv);
  gsl_rng *rng = ini_gsl_rng(args.s); // init seed, if provided

  // process files (or stdin, if none given)
  int fd = 0;
  if (!args.num_files) {
    fd = 0;
    scanFile(fd);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      fd = open(args.files[i], 0);
      scanFile(fd);
      close(fd);
    }
  }

  free_gsl_rng(rng, args.s); // save seed, if not using provided
}
