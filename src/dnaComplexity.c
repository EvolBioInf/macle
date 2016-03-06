#include "prelude.h"

#include <unistd.h>
#include <fcntl.h>

#include "args.h"
#include "gsl_rng.h"

#include "sequenceData.h"
#include "complexity.h"
#include "lempelziv.h"

void gnuplotCode(uint32_t w, uint32_t k, int n) {
  printf(//"set terminal png; "
    "set key autotitle columnheader; "
    "set ylabel \"window complexity\"; "
    "set xlabel \"window offset (w=%u, k=%u)\"; ", w, k);
  for (int i=0; i<n; i++)
    printf("%s using 1:%d with lines", i?", ''":"plot \"-\"", i+2);
  printf(";\n");
}

void printPlot(uint32_t w, uint32_t k, size_t n, Sequence *seq, double **ys) {
	printf("i ");
	for (int i=0; i<seq->numSeq; i++)
		printf("\"%s\" ", seq->headers[i]+1);
	printf("\n");
    for (size_t j=0; j*k <= n-w; j++) {
		printf("%zu ", j*k);
		for (int i=0; i<seq->numSeq; i++)
			printf("%.4f ", ys[i][j]);
		printf("\n");
	}
}

//read sequence, do stuff
void scanFile(int fd, Args *args){
  Sequence *seq = readFasta(fd);
  size_t maxlen = maxSeqLen(seq); //this implies the domain of the plot
  size_t minlen = minSeqLen(seq); //this restricts the reasonable window sizes

  //adapt window size and interval
  size_t w = args->w;
  size_t k = args->k;
  if (w==0) w=minlen;  //default window = whole (smallest) seq.
  w = MIN(w, minlen);  //biggest window = whole (smallest) seq.
  if (k==0) k=w;       //default interval = whole (smallest) seq.
  k = MIN(k, w);       //biggest interval = window size

  //array for results
  double **ys = malloc(seq->numSeq * sizeof(double*));
  for (int i=0; i<seq->numSeq; i++)
    ys[i] = calloc(maxlen, sizeof(double));

  for(int i=0; i<seq->numSeq; i++){
    char *t = getSeq(seq,i);
    size_t n = seqLen(seq,i);

    Esa *esa = getEsa(t,n);
    Fact *mlf = mlComplexity(seq, esa, w, k);
    Fact *lzf = computeLZFact(esa);

    if(args->p) {
       printf("%s \t(ML=%.4f)\n",seq->headers[i], mlf->cNor);
       printf("ML-Factors:\n");
       printFact(mlf);
       printf("LZ-Factors:\n");
       printFact(lzf);
    }

    //calculate window complexity, TODO: does this make sense?
    for (size_t j=0; j*k <= n-w; j++) {
      double l = n;
      double facs = factorsFromTo(mlf, j*k, MIN(n, j*k+w)-1);
      ys[i][j] = (facs/l - mlf->cMin)/(mlf->cMax - mlf->cMin) * ((double)w)/l;
    }

    freeFact(mlf);
    freeFact(lzf);

    freeEsa(esa);
  }

  if (!args->p) {
  if (args->g) { //print to be directly piped into gnuplot
    gnuplotCode(w,k, seq->numSeq);
    for (int i=0; i<seq->numSeq; i++) {
      printPlot(w,k,maxlen,seq,ys);
      printf("e\n");
    }
  } else { //just print resulting data
    printPlot(w,k,maxlen,seq,ys);
  }
  }

  for (int i=0; i<seq->numSeq; i++)
    free(ys[i]);
  free(ys);
  freeSequence(seq);
}

int main(int argc, char *argv[]) {
  Args args = parseArgs(argc, argv);
  gsl_rng *rng = ini_gsl_rng(args.s); //init seed, if provided

  //process files (or stdin, if none given)
  int fd = 0;
  if(!args.num_files){
    fd = 0;
    scanFile(fd, &args);
  }else{
    for(size_t i=0; i<args.num_files; i++){
      fd = open(args.files[i],0);
      scanFile(fd, &args);
      close(fd);
    }
  }

  free_gsl_rng(rng, args.s); //save seed, if not using provided
}
