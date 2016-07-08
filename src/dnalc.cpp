#include <iostream>
#include <iomanip>
using namespace std;

#include "args.h"
#include "bench.h"
#include "complexity.h"

void printIndexInfo(ComplexityData const &dat) {
  cout << "name:\t" << dat.name << endl
        << "len:\t" << dat.len << endl
        << "gc:\t" << dat.gc << endl
        << "bad:\t" << (double)dat.numbad / dat.len << endl;
  cout << "sequences:" << endl;
  for (size_t j=0; j<dat.regions.size(); j++) {
    cout << "\tindex: " << (j+1);
    cout << "\tlen: " << dat.regions[j].second;
    cout << "\tname: " << dat.labels[j] << endl;
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
  ComplexityData dat;
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
    extractData(dat, ff);

    if (args.s && !args.p) { // just dump intermediate results and quit
      saveData(dat, nullptr);
      return;
    }
  }

  for (auto &task : args.tasks) {
    if (task.idx < -1 || task.idx >= (int64_t)dat.regions.size()) {
      cerr << "ERROR: Invalid sequence index!" << endl;
      return;
    }
    size_t reglen = task.idx >= 0 ? dat.regions[task.idx].second : dat.len;
    if ((task.start!=0 || task.end!=0) && (task.start > task.end || task.start >= reglen ||
        task.end-task.start+1 > reglen - task.start)) {
      cerr << "ERROR: Invalid range!" << endl;
      return;
    }

    size_t w = args.w;
    size_t k = args.k;
    auto ys = calcComplexities(w, k, task, dat);
    if (!args.p) {
      if (args.tasks.size()==1) {
        printResults(args.w, args.w, ys, args.g);
      } else {
        cout << task.num << "\t" << ys[0].second[0] << endl;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);
  cout << fixed << setprecision(4);

  tick();
  if (args.num_files == 0)
    processFile(nullptr); //from stdin
    else
    processFile(args.files[0]);
  tock("total time");
}
