#pragma once
#include "pfasta.h"

#include <vector>
#include <string>

struct FastaSeq {
  FastaSeq(std::string const &n, std::string const &c, std::string const &s);
  std::string name;
  std::string comment;
  std::string seq;
};

/* basic sequence type representing >= 1 entry in FASTA file */
struct FastaFile {
  FastaFile();
  FastaFile(char const *file);
  std::string filename; //filename
  std::vector<FastaSeq> seqs; // the sequences
  bool failed;
};
