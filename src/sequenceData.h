/***** sequenceData.h *******************************************************
 * Description: Header file for sequence manipulation tasks.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:25 2004.
 *****************************************************************************/ 
#ifndef SEQUENCEDATA
#define SEQUENCEDATA

#include <stdio.h>
#include <stdlib.h>
#define SEQLINE 1000       /* maximal length of one line in FASTA file; hard bound */
#define SEQBUFFER 5000000  /* define the size of the sequence buffer */
#define DICSIZE 256
/* #define BORDER 'Z' */
#define BORDER '$'
/* #define BORDER '\0' */

/* basic sequence type representing >= 1 entry in FASTA file */
typedef struct sequence{
  char *seq;                /* the sequence */
  char *id;                 /* the sequence id */
  int numSeq;               /* number of sequences represented */
  int numQuery;             /* number of query sequences */
  long *borders;             /* position of last character of sequence in seq */
  char **headers;           /* FASTA header lines */
  long len;                  /* sequence length */
  int *freqTab;             /* frequency table */
  long numQueryNuc;          /* number of nucleotides in query sequence */
  long numSbjctNuc;          /* number of nucleotides in sbjct sequence */
  long numNuc;               /* number of nucleotides in sequence */
  long queryStart;
  long queryEnd;
  long effQueryNuc;          /* number of nucleotides in query that are the starting
			      * point for shustrings that are longer than expected
			      * by chance alone */
  int *sbjctId;             /* sequence id for each sbjct position */
  double queryGc;           /* GC content of query */
  double sbjctGc;           /* GC content of sbjct */
} Sequence;

char* getSeq(Sequence *seq, size_t i);
size_t seqLen(Sequence *seq, size_t i);
size_t maxSeqLen(Sequence *seq);
size_t minSeqLen(Sequence *seq);

Sequence *revcomp(Sequence *seq);
Sequence *getNextSequence(FILE *fp);

int *getDnaDictionary(int *dic);
int *getRestrictedDnaDictionary(int *dic);

Sequence *freeSequence(Sequence *seq);
Sequence *getPermanentNextSequence(FILE *fp);
void convertToAcgt(Sequence *seq);
void resetSequenceReader();
Sequence *readFasta(int fd);
Sequence **sequence2array(Sequence *seq);
void prepareSeq(Sequence *sequence);
Sequence *catSeq(Sequence *seq1, Sequence *seq2);
void randomizeSeq(Sequence *seq);
void randomizeSbjct(Sequence *seq);
Sequence *cloneSeq(Sequence *ori);
double gcContent(Sequence *seq);

#endif
