#pragma once
#include <vector>
#include <string>
#include <utility>
#include "fastafile.h"

// All information from a sequence required to calculate complexity plots
// a file stores exactly one such object with one or more regions defined
// by the fasta sequences within the file
struct ComplexityData {
  std::string name;                       // name of sequence
  size_t len;                             // length of sequence
  double gc;                              // gc content of sequence

  // for sequence regions in joined sequence:
  std::vector<std::string> labels;                //region labels
  std::vector<std::pair<size_t, size_t>> regions; //regions (start, length)

  // for global mode we need to ignore NNN... blocks:
  size_t numbad;                               // total # of bad nucleotides
  std::vector<std::pair<size_t, size_t>> bad;  // list of bad intervals (start,end)

  std::vector<size_t> fstRegionFact;              //for each region, index of first factor
  std::vector<size_t> mlf;                // match factors
};

const size_t MAX_LABEL_LEN = 32;

bool readMagic(istream &fin);
bool loadData(ComplexityData &cplx, char const *file, bool onlyInfo=false);
bool saveData(ComplexityData &cplx, char const *file);
bool renameRegions(char const *file, std::vector<std::string> const &names);

void extractData(ComplexityData &cplx, FastaFile &file);
