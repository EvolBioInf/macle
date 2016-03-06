/***** sequencedata.c *********************************************
 * Description: Collection of routines for reading and
 * manipulating sequence data.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:34:31 2004.
 *****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "sequenceData.h"
#include "stringUtil.h"
#include "eprintf.h"
#include "args.h"

static int lastSequence = 0;
static char *line = NULL;

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

/* everything that is not [acgtACGT] is flagged by a -1 */
int *getRestrictedDnaDictionary(int *dic) {
  int i;

  if (dic == NULL)
    dic = (int *)emalloc((DICSIZE + 1) * sizeof(int));

  for (i = 0; i < DICSIZE; i++)
    dic[i] = -1;

  dic['a'] = 0; /* a */
  dic['c'] = 1; /* c */
  dic['g'] = 2; /* g */
  dic['t'] = 3; /* t */
  dic['A'] = 0; /* A */
  dic['C'] = 1; /* C */
  dic['G'] = 2; /* G */
  dic['T'] = 3; /* T */

  return dic;
}

/* getDnaDictionary: create DNA dictionary */
int *getDnaDictionary(int *dic) {
  int i;

  if (dic == NULL)
    dic = (int *)malloc((DICSIZE + 1) * sizeof(int));

  for (i = 0; i < DICSIZE; i++)
    dic[i] = 0;

  dic['a'] = 0; /* a */
  dic['c'] = 1; /* c */
  dic['g'] = 2; /* g */
  dic['t'] = 3; /* t */
  dic['A'] = 0; /* A */
  dic['C'] = 1; /* C */
  dic['G'] = 2; /* G */
  dic['T'] = 3; /* T */
  dic['r'] = dic['g'];
  dic['R'] = dic['g'];
  dic['y'] = dic['t'];
  dic['Y'] = dic['t'];
  dic['m'] = dic['a'];
  dic['M'] = dic['a'];
  dic['k'] = dic['g'];
  dic['K'] = dic['g'];
  dic['s'] = dic['g'];
  dic['S'] = dic['g'];
  dic['w'] = dic['a'];
  dic['W'] = dic['a'];
  dic['h'] = dic['a'];
  dic['H'] = dic['a'];
  dic['b'] = dic['g'];
  dic['B'] = dic['g'];
  dic['v'] = dic['g'];
  dic['V'] = dic['g'];
  dic['d'] = dic['g'];
  dic['D'] = dic['g'];
  dic['n'] = dic['g'];
  dic['N'] = dic['g'];
  dic['u'] = dic['t'];
  dic['U'] = dic['t'];

  return dic;
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
  newSeq->sbjctId = NULL;
  newSeq->numSeq = 1;
  newSeq->headers = (char **)emalloc(newSeq->numSeq * sizeof(char *));
  newSeq->borders = (long *)emalloc(newSeq->numSeq * sizeof(long));
  for (i = 0; i < newSeq->numSeq; i++)
    newSeq->headers[i] = NULL;
  newSeq->id = strdup2(seq->id);
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
/* Get next sequence from an open data stream in FASTA format; this stream may be the
 * stdin */
Sequence *getPermanentNextSequence(FILE *fp) {
  Sequence *sequence;
  int seqlen, seqi, i, l;
  int currentBuffer;

  if (lastSequence) {
    return NULL;
  }
  if (line == NULL) {
    line = (char *)emalloc((SEQLINE + 2) * sizeof(char));
    line = fgets(line, SEQLINE, fp);
  }
  /* make a sequence object */
  sequence = (Sequence *)emalloc(sizeof(Sequence));
  /* allocate memory for sequence id */
  sequence->id = (char *)emalloc((strlen(line) + 1) * sizeof(char));
  /* copy sequence id */
  strcpy(sequence->id, chomp(line));
  /* allocate memory for sequence string */
  sequence->seq = (char *)emalloc((SEQBUFFER + 1) * sizeof(char));
  seqlen = 0;
  currentBuffer = SEQBUFFER;
  seqi = 0;
  while ((line = fgets(line, SEQLINE, fp)) != NULL) {
    if (strstr(line, ">") != NULL) {
      sequence->seq[seqi++] = '\0';
      sequence->seq = (char *)realloc(sequence->seq, seqi * sizeof(char));
      return sequence;
    }
    if (strlen(line) > SEQLINE) {
      printf("error in getNextSequence: cannot deal with lines longer than %d bp.\n",
             SEQLINE);
      printf("  change the SEQLINE parameter in file sequenceData.h and recompile.\n");
      exit(2);
    }
    l = strlen(line);
    /* disregard the final carriage return */
    if (line[l - 1] == '\n')
      l--;
    seqlen += l;
    if (seqlen > currentBuffer) {
      currentBuffer += SEQBUFFER;
      sequence->seq = (char *)erealloc(sequence->seq, currentBuffer);
    }
    for (i = 0; i < l; i++) {
      sequence->seq[seqi++] = line[i];
    }
    /* sequence->seq = strncat(sequence->seq,line,strlen(line)-1); */
  }
  sequence->seq[seqi++] = '\0';
  sequence->seq = (char *)realloc(sequence->seq, seqi * sizeof(char));
  sequence->len = seqi - 1;
  lastSequence = 1;
  return sequence;
}

void resetSequenceReader() {
  line = NULL;
  lastSequence = 0;
}

/* convert multiple sequences contained in seq into an
 * array of sequences each representing a single sequence
 */
Sequence **sequence2array(Sequence *seq) {
  Sequence **seqs;
  long i, j, k, len;

  /* allocate space for sequences */
  seqs = (Sequence **)emalloc(seq->numSeq * sizeof(Sequence *));
  for (i = 0; i < seq->numSeq; i++) {
    seqs[i] = (Sequence *)emalloc(sizeof(Sequence));
    seqs[i]->freqTab = (int *)emalloc(DICSIZE * sizeof(int));
    seqs[i]->seq = NULL;
    seqs[i]->borders = NULL;
    seqs[i]->numNuc = 0;
    seqs[i]->sbjctId = NULL;
    for (j = 0; j < DICSIZE; j++)
      seqs[i]->freqTab[j] = 0;
  }
  /* deal with first sequence */
  len = seq->borders[0];
  seqs[0]->seq = (char *)emalloc((len + 2) * sizeof(char));
  seqs[0]->len = len + 1;
  for (i = 0; i < len; i++) {
    seqs[0]->seq[i] = seq->seq[i];
    seqs[0]->freqTab[(int)seq->seq[i]]++;
    seqs[0]->numNuc++;
  }
  seqs[0]->numNuc *= 2;
  seqs[0]->seq[len] = BORDER;
  seqs[0]->seq[len + 1] = '\0';
  seqs[0]->id = (char *)emalloc(6 * sizeof(char));
  seqs[0]->id[0] = '\0';
  strcat(seqs[0]->id, "strId");
  seqs[0]->numSeq = 1;
  seqs[0]->borders = (long *)emalloc(sizeof(long));
  seqs[0]->borders[0] = seq->borders[0];
  seqs[0]->headers = (char **)emalloc(sizeof(char *));
  seqs[0]->headers[0] = (char *)emalloc((strlen(seq->headers[0]) + 1) * sizeof(char));
  seqs[0]->headers[0] = strcpy(seqs[0]->headers[0], seq->headers[0]);
  /* deal with remaining sequences */
  for (i = 1; i < seq->numSeq; i++) {
    len = seq->borders[i] - seq->borders[i - 1];
    seqs[i]->len = len;
    seqs[i]->seq = (char *)emalloc((len + 1) * sizeof(char));
    k = 0;
    for (j = seq->borders[i - 1] + 1; j < seq->borders[i]; j++) {
      seqs[i]->seq[k++] = seq->seq[j];
      seqs[i]->freqTab[(int)seq->seq[j]]++;
      seqs[i]->numNuc++;
    }
    seqs[i]->seq[len - 1] = BORDER;
    seqs[i]->seq[len] = '\0';
    seqs[i]->id = (char *)emalloc(6 * sizeof(char));
    seqs[i]->id[0] = '\0';
    strcat(seqs[i]->id, "strId");
    seqs[i]->numSeq = 1;
    seqs[i]->borders = (long *)emalloc(sizeof(long));
    seqs[i]->borders[0] = len - 1;
    seqs[i]->headers = (char **)emalloc(sizeof(char *));
    seqs[i]->headers[0] = (char *)emalloc((strlen(seq->headers[i]) + 1) * sizeof(char));
    seqs[i]->headers[0][0] = '\0';
    seqs[i]->headers[0] = strcpy(seqs[i]->headers[0], seq->headers[i]);
    seqs[i]->numNuc *= 2;
  }
  return seqs;
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
  s->id = (char *)emalloc(6 * sizeof(char));
  s->id[0] = '\0';
  strcat(s->id, "strId");
  headerOpen = 0;
  s->len = 0;
  s->numSeq = 0;
  s->sbjctId = NULL;
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

/* Get next sequence from an open data stream in FASTA format; this stream may be the
 * stdin */
Sequence *getNextSequence(FILE *fp) {
  Sequence *sequence;
  int seqi, i, l;
  int currentBuffer;

  if (fp == NULL || lastSequence) {
    return NULL;
  }

  if (line == NULL) {
    line = (char *)malloc((SEQLINE + 2) * sizeof(char));
    line = fgets(line, SEQLINE, fp);
  }

  /* make a sequence object */
  sequence = (Sequence *)malloc(sizeof(Sequence));
  /* allocate memory for sequence id */
  sequence->id = (char *)malloc((strlen(line) + 1) * sizeof(char));
  /* copy sequence id */
  strcpy(sequence->id, line);
  /* allocate memory for sequence string */
  sequence->seq = (char *)malloc((SEQBUFFER + 1) * sizeof(char));
  sequence->numSeq = 1;
  sequence->headers = (char **)emalloc(sizeof(char *));
  sequence->headers[0] = (char *)emalloc(sizeof(char));
  sequence->headers[0][0] = '\0';
  sequence->borders = (long *)emalloc(sizeof(long));

  sequence->len = 0;
  currentBuffer = SEQBUFFER;
  seqi = 0;
  while ((line = fgets(line, SEQLINE, fp)) != NULL) {
    if (line[0] == '>') {
      sequence->len++;
      sequence->seq = (char *)realloc(sequence->seq, sequence->len);
      sequence->seq[sequence->len - 1] = BORDER;
      sequence->borders[0] = sequence->len - 1;
      return sequence;
    }
    if (strlen(line) > SEQLINE) {
      printf("error in getNextSequence: cannot deal with lines longer than %d bp.\n",
             SEQLINE);
      printf("  change the SEQLINE parameter in file sequenceData.h and recompile.\n");
      exit(2);
    }
    l = strlen(line);
    /* disregard the final carriage return */
    if (line[l - 1] == '\n')
      l--;
    sequence->len += l;
    if (sequence->len > currentBuffer) {
      currentBuffer += SEQBUFFER;
      sequence->seq = (char *)erealloc(sequence->seq, currentBuffer);
    }
    for (i = 0; i < l; i++) {
      sequence->seq[seqi++] = line[i];
    }
    /* sequence->seq = strncat(sequence->seq,line,strlen(line)-1); */
  }
  sequence->len++;
  sequence->seq = (char *)realloc(sequence->seq, sequence->len);
  sequence->seq[sequence->len - 1] = BORDER;
  sequence->borders[0] = sequence->len - 1;
  lastSequence = 1;
  return sequence;
}

/* freeSequence: free the data structure Sequence */
Sequence *freeSequence(Sequence *seq) {
  int i;

  for (i = 0; i < seq->numSeq; i++)
    free(seq->headers[i]);
  free(seq->headers);
  free(seq->borders);
  free(seq->id);
  free(seq->seq);
  free(seq->freqTab);
  free(seq->sbjctId);
  free(seq);
  return seq;
}

/* prepareSeq: prepares sequence string for analysis by shustring-type programs.
 * Does the following: 1) set all residues to upper case
 *                     2) generate reverse complement
 *                     3) concatenate reverse complement to end of forward strand
 */
void prepareSeq(Sequence *sequence) {
  Sequence *rstrand;
  int i;

  strtoupper(sequence->seq, sequence->len);
  /* take care of reverse strand */
  rstrand = revcomp(sequence);
  rstrand->headers = (char **)emalloc(sizeof(char *));
  rstrand->headers[0] = (char *)emalloc(sizeof(char));
  rstrand->borders = (long *)emalloc(sizeof(long));
  rstrand->numSeq = 1;
  /* sequence->seq[sequence->len] = '\0'; */
  sequence->len += sequence->len;
  sequence->seq = (char *)erealloc(sequence->seq, (sequence->len + 1) * sizeof(char));
  sequence->borders =
      (long *)erealloc(sequence->borders, 2 * sequence->numSeq * sizeof(long));
  for (i = 1; i < sequence->numSeq; i++) {
    sequence->borders[2 * sequence->numSeq - i - 1] =
        sequence->len - sequence->borders[i - 1] - 2;
  }
  sequence->borders[2 * sequence->numSeq - 1] = sequence->len - 1;
  /* move first border of reverted sequences to the end */
  rstrand->seq++;
  strncat(sequence->seq, rstrand->seq, sequence->len);
  rstrand->seq--;
  sequence->seq[sequence->len - 1] = BORDER;
  sequence->seq[sequence->len] = '\0';
  freeSequence(rstrand);
}

/* catSeq: concatenate the sequences contained in two Sequence objects */
Sequence *catSeq(Sequence *seq1, Sequence *seq2) {
  Sequence *cat;
  long i, j, n;

  cat = (Sequence *)emalloc(sizeof(Sequence));
  cat->seq = (char *)emalloc((strlen(seq1->seq) + strlen(seq2->seq) + 1) * sizeof(char));
  cat->seq[0] = '\0';
  cat->seq = strncat(cat->seq, seq1->seq, seq1->len);
  cat->seq = strncat(cat->seq, seq2->seq, seq2->len);
  cat->id = (char *)emalloc(6 * sizeof(char));
  cat->id[0] = '\0';
  strcat(cat->id, "strId");
  n = seq1->numSeq + seq2->numSeq;
  cat->numSeq = n;
  cat->numQuery = seq1->numSeq;
  cat->freqTab = NULL;
  cat->queryStart = 0;
  cat->queryEnd = seq1->len - 1;
  cat->borders = (long *)emalloc(2 * n * sizeof(long));
  cat->headers = (char **)emalloc(n * sizeof(char *));
  /* take care of the n headers */
  for (i = 0; i < seq1->numSeq; i++) {
    cat->headers[i] = (char *)emalloc((strlen(seq1->headers[i]) + 1) * sizeof(char));
    cat->headers[i] = strcpy(cat->headers[i], seq1->headers[i]);
  }
  j = i;
  for (i = 0; i < seq2->numSeq; i++) {
    cat->headers[j] = (char *)emalloc((strlen(seq2->headers[i]) + 1) * sizeof(char));
    cat->headers[j] = strcpy(cat->headers[j], seq2->headers[i]);
    j++;
  }
  /* take care of the 2n borders */
  for (i = 0; i < seq1->numSeq; i++) {
    cat->borders[i] = seq1->borders[i];
  }
  j = i;
  for (i = 0; i < seq2->numSeq; i++) {
    cat->borders[j + i] = seq1->borders[2 * seq1->numSeq - 1] + seq2->borders[i] + 1;
  }
  /* sbjct IDs */
  cat->sbjctId = (int *)emalloc(seq2->len * sizeof(int));
  for (i = 0; i < seq2->len; i++)
    cat->sbjctId[i] = -1;
  for (i = 0; i <= seq2->borders[0]; i++)
    cat->sbjctId[i] = 0;
  for (i = 1; i < seq2->numSeq; i++)
    for (j = seq2->borders[i - 1] + 1; j <= seq2->borders[i]; j++)
      cat->sbjctId[j] = i;
  n = seq2->numSeq - 1;
  for (; i < 2 * seq2->numSeq; i++) {
    for (j = seq2->borders[i - 1] + 1; j <= seq2->borders[i]; j++)
      cat->sbjctId[j] = n;
    n--;
  }

  cat->len = seq1->len + seq2->len;
  cat->numQueryNuc = seq1->numNuc;
  cat->numSbjctNuc = seq2->numNuc;
  cat->numNuc = seq1->numNuc + seq2->numNuc;
  return cat;
}

/* cloneSeq: make exact copy of Sequence object */
Sequence *cloneSeq(Sequence *ori) {
  Sequence *clone;
  int i;

  clone = (Sequence *)emalloc(sizeof(Sequence));
  clone->sbjctId = NULL;
  clone->seq = (char *)emalloc(((int)ori->len) * sizeof(char));
  clone->seq[0] = '\0';
  clone->seq = strncpy(clone->seq, ori->seq, ori->len);
  /* clone->seq[ori->len] = '\0'; */
  clone->id = (char *)emalloc(6 * sizeof(char));
  clone->id = strncpy(clone->id, ori->id, 6);
  clone->numSeq = ori->numSeq;
  clone->numQuery = ori->numQuery;
  clone->borders = (long *)emalloc(ori->numSeq * sizeof(long));
  for (i = 0; i < ori->numSeq; i++)
    clone->borders[i] = ori->borders[i];
  clone->headers = (char **)emalloc(ori->numSeq * sizeof(char *));
  for (i = 0; i < ori->numSeq; i++) {
    clone->headers[i] = (char *)emalloc((strlen(ori->headers[i]) + 1) * sizeof(char));
    clone->headers[i] = strcpy(clone->headers[i], ori->headers[i]);
  }
  clone->len = ori->len;
  clone->freqTab = (int *)emalloc(DICSIZE * sizeof(int));
  for (i = 0; i < DICSIZE; i++)
    clone->freqTab[i] = ori->freqTab[i];
  clone->numQueryNuc = ori->numQueryNuc;
  clone->numSbjctNuc = ori->numSbjctNuc;
  clone->numNuc = ori->numNuc;
  clone->queryStart = ori->queryStart;
  clone->queryEnd = ori->queryEnd;

  return clone;
}

// helpers

double gcContent(Sequence *seq) {
  uint64_t gc = 0;
  uint64_t numChar = 0;
  for (int i = 0; i < seq->numSeq; i++) {
    char *s = getSeq(seq, i);
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

char *getSeq(Sequence *seq, size_t i) {
  return &(seq->seq[i ? seq->borders[i - 1] + 2 : 0]);
}

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
size_t minSeqLen(Sequence *seq) {
  size_t min = seq->len;
  for (int i = 0; i < seq->numSeq; i++) {
    size_t len = seqLen(seq, i);
    if (len < min)
      min = len;
  }
  return min;
}
