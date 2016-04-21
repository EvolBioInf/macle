#include "prelude.h"

#include "args.h"
#include "bench.h"
#include "eprintf.h"
#include "gsl_rng.h"

#include "fastafile.h"
#include "matchlength.h"
#include "lempelziv.h"
#include "list.h"
#include "periodicity.h"
#include "complexity.h"

// calculate the GC content
double gcContent(const pfasta_seq *ps) {
  size_t gc = 0;
  char const *ptr = ps->seq;
  for (; *ptr; ptr++)
    if (*ptr == 'g' || *ptr == 'G' || *ptr == 'c' || *ptr == 'C')
      gc++;
  return (double)gc / (ptr - ps->seq);
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

void printPlot(uint32_t k, size_t n, FastaFile *ff, double **ys) {
  printf("offset ");
  for (size_t i = 0; i < ff->n; i++) { // two columns for each seq
    printf("\"%s %s\" ", ff->seq[i].name, "(ML)");
    printf("\"%s %s\" ", ff->seq[i].name, "(Per)");
  }
  printf("\n");
  for (size_t j = 0; j < n; j++) {
    printf("%zu ", j * k);
    for (size_t i = 0; i < 2 * ff->n; i++)
      printf("%.4f ", ys[i][j]);
    printf("\n");
  }
}

// read sequence, do stuff
void scanFile(FastaFile *ff) {
  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  for (size_t i = 0; i < ff->n; i++) {
    size_t len = ff->seq[i].len;
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
  double **ys = (double**)emalloc(2 * ff->n * sizeof(double *));
  for (size_t i = 0; i < 2 * ff->n; i++)
    ys[i] = (double*)ecalloc(entries + 1, sizeof(double));

  for (size_t i = 0; i < ff->n; i++) {
    double gc = gcContent(&(ff->seq[i]));
    char *t = ff->seq[i].seq;
    size_t n = ff->seq[i].len;

    tick();
    Esa *esa = getEsa(t, n + 1); // esa for sequence+$
    tock("getEsa");
    if (args.p) {
      printf("%s %s\n", ff->seq[i].name, ff->seq[i].comment);
      /* printEsa(esa); */
    }

    tick();
    Fact *mlf = computeMLFact(esa);
    tock("computeMLFact");
    if (args.p) {
      printf("ML-Factors (%zu):\n", mlf->n);
      printFact(mlf);
    }

    tick();
    mlComplexity(w, k, ys[2 * i], mlf, gc);
    tock("mlComplexity");
    freeFact(mlf);

    tick();
    Fact *lzf = computeLZFact(esa, false);
    tock("computeLZFact");
    tick();
    size_t pnum;
    List **ls = getPeriodicityLists(true, lzf, esa, &pnum);
    tock("getPeriodicityLists");

    tick();
    runComplexity(w, k, ys[2 * i + 1], n + 1, ls);
    tock("runComplexity");

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
      if (args.gf==0) {
        for (size_t i = 0; i < 2 * ff->n; i++) {
          for (size_t j = 0; j < entries; j++)
            printf("%zu %.4f\n", j * k, ys[i][j]);
          printf("\n");
        }
      } else if (args.gf==1) {
        gnuplotCode(w, k, 2 * ff->n);
        for (size_t i = 0; i < 2 * ff->n; i++) {
          printPlot(k, entries, ff, ys);
          printf("e\n");
        }
      }
    } else { // just print resulting data
      printPlot(k, entries, ff, ys);
    }
  }

  for (size_t i = 0; i < 2 * ff->n; i++)
    free(ys[i]);
  free(ys);
}

int main(int argc, char *argv[]) {
  progname = PROGNAME;
  parseArgs(argc, argv);
  gsl_rng *rng = ini_gsl_rng(args.s); // init seed, if provided

  // process files (or stdin, if none given)
  if (!args.num_files) {
    FastaFile *ff = read_fasta_file(NULL);
    scanFile(ff);
    free_fasta_file(ff);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      tick();
      FastaFile *ff = read_fasta_file(args.files[i]);
      tock("readFastaFromFile");
      if (!ff) {
        fprintf(stderr, "Skipping invalid FASTA file...\n");
        continue;
      }
      scanFile(ff);
      free_fasta_file(ff);
    }
  }

  free_gsl_rng(rng, args.s); // save seed, if not using provided
}
