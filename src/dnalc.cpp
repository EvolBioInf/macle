#include <iostream>
#include <algorithm>
#include <iomanip>
#include <queue>
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
    cout << (i ? ", ''" : "plot \"$PLOTFILE\"") << " using 2:"<<(i+3)<<" with lines";
  cout << ";" << endl;
}

// print data: X Y1 ... Yn
void printPlot(int64_t idx, vector<pair<size_t,size_t>> &regs, uint32_t w, uint32_t k, ResultMat const &ys) {
  queue<pair<size_t,size_t>> rs;
  for (auto &r : regs)
    rs.emplace(r);
  uint32_t rcnt = idx < 0 ? 1 : idx+1; //region counter for output

  for (size_t j = 0; j < ys[0].second.size(); j++) {
    size_t off = j * k + w / 2;
    // cout << regs[idx].first + off << "\t";
    if (idx < 0) {
      if (off >= rs.front().first + rs.front().second) {
        rs.pop();
        rcnt++;
      }
      off -= rs.front().first;
    }
    cout << rcnt << "\t" << off << "\t"; // center of window
    for (size_t i = 0; i < ys.size(); i++)
      cout << ys[i].second[j] <<"\t";
    cout << endl;
  }
}

void printResults(int64_t idx, vector<pair<size_t,size_t>> &regs, size_t w, size_t k, ResultMat const &ys, bool gnuplot) {
  if (!gnuplot) { //simple output
    printPlot(idx, regs, w, k, ys);
  } else { //dnalc_plot
    cout << "DNALC_PLOT" << endl; //magic keyword
    gnuplotCode(w, k, ys.size()); // gnuplot control code
    // print column header (for plot labels)
    cout << "offset\t";
    for (size_t j = 0; j < ys.size(); j++) // columns for each seq
      cout << "\"" << ys[j].first << "\"\t";
    cout << endl;
    printPlot(idx, regs, w, k, ys); // print plot itself
  }
}

bool check_unique_names(FastaFile &ff) {
  vector<string> names;
  for (auto &seq : ff.seqs)
    names.push_back(seq.name.substr(0,MAX_LABEL_LEN));
  sort(names.begin(),names.end());
  for (size_t i=0; i<names.size()-1; i++)
    if (names[i]==names[i+1])
      return false;
  return true;
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
      cerr << "Invalid FASTA file!" << endl;
      return;
    }
    if (!check_unique_names(ff)) {
      cerr << "Headers of the FASTA sequence must be unique before the first whitespace or 32 characters!" << endl;
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
        printResults(task.idx, dat.regions, w, k, ys, args.g);
      } else {
        cout << task.num << "\t" << ys[0].second[0] << endl;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);
  cout << fixed << setprecision(4);

  if (args.newnames.size()>0) {
    if (args.num_files==0) {
      cerr << "ERROR: No index file provided!" << endl;
      return EXIT_FAILURE;
    }
    renameRegions(args.files[0], args.newnames);
    return EXIT_SUCCESS;
  }

  tick();
  if (args.num_files == 0)
    processFile(nullptr); //from stdin
  else
    processFile(args.files[0]);
  tock("total time");
}
