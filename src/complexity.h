#pragma once
#include <vector>
#include <list>
#include <utility>

#include "args.h"
#include "index.h"

size_t numEntries(size_t n, size_t w, size_t k);

void mlComplexity(size_t offset, size_t n, size_t w, size_t k, std::vector<double> &y, ComplexityData const &dat);

typedef std::vector<std::pair<std::string,std::vector<double>>> ResultMat;

// This complicated function calculates the complexity data depending on mode.
// The input data can be "joined" -> one single sequence with "regions", or all
// sequences in the input are separate.
// Also the user can use global mode (one value per complexity and sequence),
// and the user can choose a region or sequence to process.
// If NOT joined: no chosen seqnum -> compute for all separately, otherwise only given sequence
// If joined: no chosen seqnum -> compute for complete sequence, otherwise only one region
ResultMat calcComplexities(size_t &w, size_t &k, Task task, ComplexityData const &dat);
