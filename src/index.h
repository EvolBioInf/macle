#include <vector>
#include <string>
#include <utility>
#include "periodicity.h"

// All information from a sequence required to calculate complexity plots
struct ComplexityData {
  std::string name;                       // name of sequence
  size_t len;                             // length of sequence
  double gc;                              // gc content of sequence
  std::vector<size_t> mlf;                // match factors
  std::vector<std::list<Periodicity>> pl; // periodicities
  std::vector<std::pair<size_t, size_t>> bad;  // list of bad intervals
};

bool loadData(std::vector<ComplexityData> &cplx, char const *file);
void saveData(std::vector<ComplexityData> &vec);
