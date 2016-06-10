#include <iostream>
#include <iomanip>
using namespace std;

#include "args.h"
#include "bench.h"
#include "complexity.h"

void printIndexInfo(vector<ComplexityData> const &dat) {
  for (size_t i=0; i<dat.size(); i++) {
    if (dat.size()>1)
      cout << (i+1) << ":" << endl;
    string tab = dat.size()>1 ? "\t" : "";
    cout << tab << "name:\t" << dat[i].name << endl
          << tab << "len:\t" << dat[i].len << endl
          << tab << "gc:\t" << dat[i].gc << endl
          << tab << "bad:\t" << (double)dat[i].numbad / dat[i].len << endl;
    if (dat[i].regions.size()>1) {
      cout << tab << "regions:" << endl;
      for (size_t j=0; j<dat[i].regions.size(); j++) {
        auto &r = dat[i].regions[j];
        cout << "\tindex:\t" << (j+1) << ":\tregion:\t" << r.first << "-" << r.first+r.second-1 << endl;
        cout << "\tname:\t" << dat[i].labels[j] << endl;
      }
    }
  }
}

void gnuplotCode(uint32_t w, uint32_t k, int n) {
  cout << "set key autotitle columnheader; set ylabel \"window complexity\"; "
    << "set xlabel \"window offset (w=" << w << ", k="<<k<<")\"; ";
  for (int i = 0; i < n; i++)
    cout << (i ? ", ''" : "plot \"$PLOTFILE\"") << " using 1:"<<(i+2)<<" with lines";
  cout << ";" << endl;
}

// print data: X Y1 ... Yn
void printPlot(uint32_t w, uint32_t k, ResultMat const &ys) {
  for (size_t j = 0; j < ys[0].second.size(); j++) {
    cout << (j * k + w / 2) << "\t"; // center of window
    for (size_t i = 0; i < ys.size(); i++)
      cout << ys[i].second[j] <<"\t";
    cout << endl;
  }
}

void printResults(size_t w, size_t k, ResultMat const &ys, int mode) {
  cout << fixed << setprecision(4);
  if (mode == 0) { //simple output
    printPlot(w, k, ys);
  } else if (mode == 2) {
    for (size_t i = 0; i < ys.size(); i++) {
      for (size_t j = 0; j < ys[0].second.size(); j++)
        cout << (j * k + w / 2) << "\t" << ys[i].second[j] << endl;
      cout << endl;
    }
  } else if (mode == 1) { //dnalc_plot
    cout << "DNALC_PLOT" << endl; //magic keyword
    gnuplotCode(w, k, ys.size()); // gnuplot control code
    // print column header (for plot labels)
    cout << "offset\t";
    for (size_t j = 0; j < ys.size(); j++) // columns for each seq
      cout << "\"" << ys[j].first << "\"\t";
    cout << endl;
    printPlot(w, k, ys); // print plot itself
  }
}

//load / extract data, show results
void processFile(char const *file) {
  vector<ComplexityData> dat;
  if (args.i) { //load from index
    tick();
    if (!loadData(dat, file, args.l))
      return;
    tock("loadData");
    if (args.l) { //list index file contents and exit
      printIndexInfo(dat);
      return;
    }
  } else { // not loading from pre-computed data -> fasta file
    tick();
    FastaFile ff(file);
    tock("readFastaFromFile");
    if (ff.failed) {
      cerr << "Skipping invalid FASTA file..." << endl;
      return;
    }
    extractData(dat, ff, args.j);

    if (args.s && !args.p) { // just dump intermediate results and quit
      saveData(dat, nullptr);
      return;
    }
  }

  if (args.n>max(dat.size(),dat[0].regions.size())) {
    cerr << "Invalid sequence number!" << endl;
    return;
  }

  size_t w = args.w;
  size_t k = args.k;
  auto ys = calcComplexities(w, k, args.m, args.n, dat);
  if (!args.p)
    printResults(w, k, ys, args.g);
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  tick();
  if (!args.num_files) {
    processFile(nullptr); //from stdin
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      processFile(args.files[i]);
    }
  }
  tock("total time");
}
