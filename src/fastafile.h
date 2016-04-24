#pragma once
#include "prelude.h"
#include "pfasta.h"

#include <vector>

/* basic sequence type representing >= 1 entry in FASTA file */
class FastaFile {
public:
  FastaFile(char const *file);
  ~FastaFile();

  std::vector<pfasta_seq> seqs; /* the sequences */
  bool failed;
};
