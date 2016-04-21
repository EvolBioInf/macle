#pragma once
#include "prelude.h"
#include "pfasta.h"

/* basic sequence type representing >= 1 entry in FASTA file */
typedef struct fastafile {
  pfasta_seq *seq; /* the sequences */
  size_t n;        /* number of sequences */
} FastaFile;

FastaFile *read_fasta_file(char const *file);
void free_fasta_file(FastaFile *fastafile);
