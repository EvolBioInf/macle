#include <iostream>
#include <algorithm>
#include <iomanip>
#include <queue>
#include <map>
using namespace std;

#include "args.h"
#include "bench.h"
#include "complexity.h"
#include "util.h"

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
void printPlot(Task &t, vector<string> &lbls, vector<pair<size_t,size_t>> &regs, uint32_t w, uint32_t k, ResultMat const &ys) {
  queue<pair<size_t,size_t>> rs;
  for (auto &r : regs)
    rs.emplace(r);
  uint32_t rcnt = t.idx < 0 ? 0 : t.idx; //region counter for output

  for (size_t j = 0; j < ys[0].second.size(); j++) {
    size_t off = j * k + w / 2;
    // cout << regs[idx].first + off << "\t";
    if (t.idx < 0) {
      if (off >= rs.front().first + rs.front().second) {
        rs.pop();
        rcnt++;
      }
      off -= rs.front().first;
    }
    string lbl = (t.idx<0 && ys[0].second.size()==1) ? "<file>" : lbls[rcnt];
    cout << lbl << "\t" << off << "\t"; // center of window
    for (size_t i = 0; i < ys.size(); i++)
      cout << ys[i].second[j] <<"\t";
    cout << endl;
  }
}

void printResults(Task &t, vector<string> &lbls, vector<pair<size_t,size_t>> &regs, size_t w, size_t k, ResultMat const &ys, bool gnuplot) {
  if (!gnuplot) { //simple output
    printPlot(t, lbls, regs, w, k, ys);
  } else { //macle_plot
    cout << "MACLE_PLOT" << endl; //magic keyword
    gnuplotCode(w, k, ys.size()); // gnuplot control code
    // print column header (for plot labels)
    cout << "offset\t";
    for (size_t j = 0; j < ys.size(); j++) // columns for each seq
      cout << "\"" << ys[j].first << "\"\t";
    cout << endl;
    printPlot(t, lbls, regs, w, k, ys); // print plot itself
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

  //infer whether given file is an index (user can forget -i)
  if (file && with_file_in(file, readMagic))
    args.i = true;

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

  //map from region name to index within index file
  map<string, int64_t> nameidx;
  nameidx[""] = -1;
  for (int64_t i=0; i<(int64_t)dat.labels.size(); i++)
    nameidx[dat.labels[i]] = i;

  for (auto &task : args.tasks) {
    string errstr = "ERROR in task #" + to_string(task.num) + ": ";

    //get index of region if not global adressing
    if (task.lbl != "") {
      if (nameidx.find(task.lbl) == nameidx.end()) {
        cerr << errstr << "Invalid sequence name: " << task.lbl << endl;
        continue;
      }
      task.idx = nameidx[task.lbl];
    }
    //sanity check for manually set index values
    if (task.idx < -1 || task.idx >= (int64_t)dat.regions.size()) {
      cerr << errstr << "Invalid sequence index (" << task.idx << ") for name " << task.lbl << endl;
      continue;
    }

    size_t reglen = task.idx >= 0 ? dat.regions[task.idx].second : dat.len;
    if ((task.start!=0 || task.end!=0) && (task.start > task.end || task.start >= reglen ||
        task.end-task.start+1 > reglen - task.start)) {
      cerr << errstr << "Invalid range: " << task.start << "-" << task.end << endl;
      continue;
    }

    size_t w = args.w;
    size_t k = args.k;
    auto ys = calcComplexities(w, k, task, dat);
    if (!args.p) {
      if (args.tasks.size()==1) {
        printResults(task, dat.labels, dat.regions, w, k, ys, args.g);
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
