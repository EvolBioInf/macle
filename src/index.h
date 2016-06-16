#pragma once
#include <vector>
#include <string>
#include <utility>
#include "fastafile.h"

// All information from a sequence required to calculate complexity plots
// either a file stores multiple such objects for separate sequences with one
// region, or a file stores exactly one such object with one or more regions
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

bool loadData(std::vector<ComplexityData> &cplx, char const *file, bool onlyInfo=false, size_t idx=0);
bool saveData(std::vector<ComplexityData> &cplx, char const *file);

void extractData(std::vector<ComplexityData> &cplx, FastaFile &file, bool joinSeqs);
