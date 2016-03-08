/***** sequencedata.c *********************************************
 * Description: Collection of routines for reading and
 * manipulating sequence data.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:31 2004.
 *****************************************************************/
#include "prelude.h"
#include <string.h>
#include <ctype.h>
#include <fcntl.h>
#include <unistd.h>

#include "sequenceData.h"
#include "eprintf.h"
#include "stringUtil.h"

/* convertToACGT: convert nucleotide data to acgt alphabet.
 */
void convertToAcgt(Sequence *seq) {
  replace((char *)seq->seq, 'r', 'g');
  replace((char *)seq->seq, 'y', 't');
  replace((char *)seq->seq, 'm', 'a');
  replace((char *)seq->seq, 'k', 'g');
  replace((char *)seq->seq, 's', 'g');
  replace((char *)seq->seq, 'w', 'a');
  replace((char *)seq->seq, 'h', 'a');
  replace((char *)seq->seq, 'b', 'g');
  replace((char *)seq->seq, 'v', 'g');
  replace((char *)seq->seq, 'd', 'g');
  replace((char *)seq->seq, 'n', 'g');
  replace((char *)seq->seq, 'u', 't');

  replace((char *)seq->seq, 'R', 'G');
  replace((char *)seq->seq, 'Y', 'T');
  replace((char *)seq->seq, 'M', 'A');
  replace((char *)seq->seq, 'K', 'G');
  replace((char *)seq->seq, 'S', 'G');
  replace((char *)seq->seq, 'W', 'A');
  replace((char *)seq->seq, 'H', 'A');
  replace((char *)seq->seq, 'B', 'G');
  replace((char *)seq->seq, 'V', 'G');
  replace((char *)seq->seq, 'D', 'G');
  replace((char *)seq->seq, 'N', 'G');
  replace((char *)seq->seq, 'U', 'T');
}

/* reverse and complement a sequence */
Sequence *revcomp(Sequence *seq) {
  long i, j, n;
  char c;
  Sequence *newSeq;
  newSeq = (Sequence *)emalloc(sizeof(Sequence));

  /*   n = strlen(seq->seq); */
  n = seq->len;
  newSeq->seq = (char *)emalloc((n + 1) * sizeof(char));
  newSeq->freqTab = NULL;
  newSeq->numSeq = 1;
  newSeq->headers = (char **)emalloc(newSeq->numSeq * sizeof(char *));
  newSeq->borders = (long *)emalloc(newSeq->numSeq * sizeof(long));
  for (i = 0; i < newSeq->numSeq; i++)
    newSeq->headers[i] = NULL;
  j = 0;
  for (i = n - 1; i >= 0; i--) {
    c = seq->seq[i];
    switch (c) {
    case BORDER:
      newSeq->seq[j++] = BORDER;
      break;
    case 'A':
      newSeq->seq[j++] = 'T';
      break;
    case 'C':
      newSeq->seq[j++] = 'G';
      break;
    case 'G':
      newSeq->seq[j++] = 'C';
      break;
    case 'T':
      newSeq->seq[j++] = 'A';
      break;
    default:
      newSeq->seq[j++] = c;
      break;
    }
  }
  newSeq->seq[n] = '\0';
  return newSeq;
}

/* read FASTA-formatted sequence data from an open file descriptor
 * into single sequence string
 */
Sequence *readFasta(int fd) {
  Sequence *s;
  char buf[BUFSIZ];
  int headerOpen;
  int headerLen;
  long i, maxLen;
  int c;

  if (fd < 0)
    return NULL;

  s = (Sequence *)emalloc(sizeof(Sequence));
  s->freqTab = (int *)emalloc(DICSIZE * sizeof(int));
  for (i = 0; i < DICSIZE; i++)
    s->freqTab[i] = 0;
  s->borders = (long *)emalloc(sizeof(long));
  s->headers = (char **)emalloc(sizeof(char *));
  headerOpen = 0;
  s->len = 0;
  s->numSeq = 0;
  maxLen = 0;
  headerLen = 0;

  while ((c = read(fd, buf, BUFSIZ)) > 0) {
    if (s->len + c + 1 > maxLen) {
      if (maxLen >= BUFSIZ) {
        maxLen *= 2;
        s->seq = (char *)erealloc(s->seq, (maxLen + 2) * sizeof(char));
      } else {
        maxLen = BUFSIZ;
        s->seq = (char *)emalloc((maxLen + 2) * sizeof(char));
      }
    }
    /* go through the c characters just read into buf */
    for (i = 0; i < c; i++) {
      if (buf[i] == '>') {
        headerOpen = 1;
        /* take care of sequence borders */
        s->borders = (long *)erealloc(s->borders, (s->numSeq + 1) * sizeof(long));
        if (s->numSeq >= 1) {
          s->seq[s->len] = BORDER; /* unique character to terminate each sequence */
          s->borders[s->numSeq - 1] = s->len;
          s->len++;
          // NEW: 0-terminate in-between
          s->seq[s->len++] = '\0';
        }
        /* take care of sequence headers */
        s->headers = (char **)erealloc(s->headers, (s->numSeq + 1) * sizeof(char *));
        s->headers[s->numSeq] = (char *)emalloc((BUFSIZ + 1) * sizeof(char));
        headerLen = 0;
        s->numSeq++;
      }
      if (headerOpen) {
        if (buf[i] == '\n') {
          headerOpen = 0;
          s->headers[s->numSeq - 1][headerLen] = '\0';
          /* trim header to actual length */
          s->headers[s->numSeq - 1] =
              (char *)erealloc(s->headers[s->numSeq - 1], (headerLen + 1) * sizeof(char));
        } else if (headerLen < BUFSIZ && isprint(buf[i]))
          s->headers[s->numSeq - 1][headerLen++] = buf[i];
      } else if (buf[i] != '\n') {
        s->seq[s->len] = buf[i];
        s->freqTab[(int)buf[i]]++;
        s->len++;
      }
    }
  }
  /* add last border */
  if (s->len < maxLen)
    s->seq[s->len] = BORDER;
  else {
    printf("ERROR [readFasta]: s->len: %d; maxLen: %d\n", (int)s->len, (int)maxLen);
    exit(0);
  }
  s->len++;
  /* set end of last sequence read */
  s->borders[s->numSeq - 1] = s->len - 1;
  /* trim sequence to actual size */
  s->seq = (char *)erealloc(s->seq, (s->len + 1) * sizeof(char));
  s->seq[s->len] = '\0';
  return s;
}

/* freeSequence: free the data structure Sequence */
void freeSequence(Sequence *seq) {
  for (int i = 0; i < seq->numSeq; i++)
    free(seq->headers[i]);
  free(seq->headers);
  free(seq->borders);
  free(seq->seq);
  free(seq->freqTab);
  free(seq);
}

// helpers

// TODO: calculate gc content from freqTab directly?
double gcContent(Sequence *seq) {
  uint64_t gc = 0;
  uint64_t numChar = 0;
  for (int i = 0; i < seq->numSeq; i++) {
    char *s = seqStr(seq, i);
    size_t len = seqLen(seq, i);
    for (size_t j = 0; j < len; j++) {
      if (s[j] == 'A' || s[j] == 'C' || s[j] == 'G' || s[j] == 'T') {
        numChar++;
        if (s[j] == 'G' || s[j] == 'C')
          gc++;
      }
    }
  }
  return (double)gc / (double)(numChar);
}

// get pointer to start of i'th sequence
char *seqStr(Sequence *seq, size_t i) {
  return &(seq->seq[i ? seq->borders[i - 1] + 2 : 0]);
}
// get real length of i'th sequence (excluding trailing $)
size_t seqLen(Sequence *seq, size_t i) {
  return seq->borders[i] - (i ? seq->borders[i - 1] + 2 : 0);
}

// length of longest sequence in file
size_t maxSeqLen(Sequence *seq) {
  size_t max = 0;
  for (int i = 0; i < seq->numSeq; i++) {
    size_t len = seqLen(seq, i);
    if (len > max)
      max = len;
  }
  return max;
}
// length of shortest sequence in file
size_t minSeqLen(Sequence *seq) {
  size_t min = seq->len;
  for (int i = 0; i < seq->numSeq; i++) {
    size_t len = seqLen(seq, i);
    if (len < min)
      min = len;
  }
  return min;
}

// read sequence(s) from given file, if NULL, uses stdin
Sequence *readFastaFromFile(char *file) {
  int fd = file ? open(file, 0) : 0;
  Sequence *seq = readFasta(fd);
  if (file)
    close(fd);
  return seq;
}
