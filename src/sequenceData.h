/***** sequenceData.h *******************************************************
 * Description: Header file for sequence manipulation tasks.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:25 2004.
 *****************************************************************************/
#ifndef SEQUENCEDATA
#define SEQUENCEDATA

#include <stdio.h>
#include <stdlib.h>
#define SEQLINE 1000      /* maximal length of one line in FASTA file; hard bound */
#define SEQBUFFER 5000000 /* define the size of the sequence buffer */
#define DICSIZE 256
/* #define BORDER 'Z' */
#define BORDER '$'
/* #define BORDER '\0' */

/* basic sequence type representing >= 1 entry in FASTA file */
typedef struct sequence {
  char *seq;        /* the sequence */
  long len;         /* sequence length */

  int numSeq;       /* number of sequences represented */
  long *borders;    /* position of last character of sequence in seq */
  char **headers;   /* FASTA header lines */
  int *freqTab;     /* frequency table */
} Sequence;

Sequence *readFasta(int fd);
void freeSequence(Sequence *seq);

char *seqStr(Sequence *seq, size_t i);
size_t seqLen(Sequence *seq, size_t i);
size_t maxSeqLen(Sequence *seq);
size_t minSeqLen(Sequence *seq);

Sequence *revcomp(Sequence *seq);
void convertToAcgt(Sequence *seq);
double gcContent(Sequence *seq);

#endif
