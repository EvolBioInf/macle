#pragma once
#include "pfasta.h"

#include <vector>
#include <string>

struct FastaSeq {
  FastaSeq(std::string n, std::string c, std::string s);
  std::string name;
  std::string comment;
  std::string seq;
};

/* basic sequence type representing >= 1 entry in FASTA file */
struct FastaFile {
  FastaFile(char const *file);
  std::vector<FastaSeq> seqs; /* the sequences */
  bool failed;
};
