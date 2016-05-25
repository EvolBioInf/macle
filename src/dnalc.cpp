#include <iostream>
#include <iomanip>
#include <utility>
using namespace std;

#include "args.h"
#include "bench.h"

#include "fastafile.h"
#include "index.h"
#include "matchlength.h"
#include "lempelziv.h"
#include "periodicity.h"
#include "complexity.h"
#include "util.h"

using namespace std;

// given sequences from a fasta file, calculate match factors and runs
void extractData(vector<ComplexityData> &cplx, FastaFile &file, bool joinSeqs) {
  FastaFile *ff = &file;

  FastaFile virtfile;
  if (joinSeqs) { //construct concatenated sequence
    string wholeseq = "";
    for (auto &seq : file.seqs) {
      wholeseq += seq.seq;
    }
    virtfile.seqs.push_back(FastaSeq(file.filename, "", wholeseq));
    ff = &virtfile;
  }

  for (auto &seq : ff->seqs) {
    ComplexityData c;
    c.name = seq.name;
    c.gc = gcContent(seq.seq);
    c.len = seq.seq.size();

    string &s = seq.seq;
    s = s + "$" + revComp(s) + "$";
    tick();
    Esa esa(s.c_str(), s.size()); // esa for seq+$+revseq+$
    tock("getEsa (2n)");
    tick();
    Fact mlf;
    computeMLFact(mlf, esa);
    tock("computeMLFact");
    c.mlf = mlf.fact;

    if (args.p) {
      cout << seq.name << seq.comment << endl;
      // esa.print();
    }

    if (args.p) {
      cout << "ML-Factors for both strands (" << mlf.fact.size() << "):" << endl;
      mlf.print();
    }

    tick();
    s.resize(s.size() / 2); // drop complementary seq.
    s.shrink_to_fit();
    esa.str = s.c_str();
    esa.n = s.size();
    reduceEsa(esa);
    tock("reduceEsa");


    tick();
    Fact lzf;
    computeLZFact(lzf, esa, false);
    lzf.lpf.resize(0); // we dont need it, memory waste
    tock("computeLZFact");
    tick();
    size_t pnum;
    c.pl = getPeriodicityLists(true, lzf, esa, pnum);
    tock("getPeriodicityLists");

    if (args.p) {
      cout << "LZ-Factors (" << lzf.fact.size() << "):" << endl;
      lzf.print();
      cout << "Runs (" << pnum << "):" << endl;
      for (auto &l : c.pl)
        for (auto p : l)
          printPeriodicity(p);
    }

    if (joinSeqs) { //store region information in concat. seq
      size_t offset=0;
      for (auto &it : file.seqs) {
        c.regions.push_back(make_pair(offset, it.seq.size()));
        c.labels.push_back(it.name+" "+it.comment);
        offset += it.seq.size();
      }
    }

    // get list of bad intervals
    tick();
    bool insidebad = false;
    size_t start = 0;
    string validchars = "ACGT$";
    for (size_t i = 0; i < s.size(); i++) {
      bool valid = validchars.find(s[i]) != string::npos;
      if (insidebad && valid) {
        insidebad = false;
        c.bad.push_back(make_pair(start, i - 1));
      } else if (!insidebad && !valid) {
        start = i;
        insidebad = true;
      }
    }
    if (insidebad) // push last one, if we are inside
      c.bad.push_back(make_pair(start, s.size() - 1));

    //calculate total number of bad nucleotides for global mode
    c.numbad=0;
    for (auto &bad : c.bad)
      c.numbad += bad.second - bad.first + 1;

    tock("find bad intervals");

    cplx.push_back(c);
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
void printPlot(uint32_t w, uint32_t k, size_t n, vector<vector<double>> &ys) {
  for (size_t j = 0; j < n; j++) {
    cout << (j * k + w / 2) << "\t"; // center of window
    for (size_t i = 0; i < ys.size(); i++)
      cout << ys[i][j] <<"\t";
    cout << endl;
  }
}

// read data, do stuff
void processFile(char const *file) {
  vector<ComplexityData> dat;
  if (args.i) { //load from index
    tick();
    if (!loadData(dat, file))
      return;
    tock("loadData");
    if (args.l) { //list index file contents and exit
      for (size_t i=0; i<dat.size(); i++) {
        if (dat.size()>1)
          cout << i << ":" << endl;
        string tab = dat.size()>1 ? "\t" : "";
        cout << tab << "name:\t" << dat[i].name << endl
             << tab << "len:\t" << dat[i].len << endl
             << tab << "gc:\t" << dat[i].gc << endl
             << tab << "bad:\t" << (double)dat[i].numbad / dat[i].len << endl;
        if (dat[i].regions.size()>1) {
          cout << tab << "regions:" << endl;
          for (size_t j=0; j<dat[i].regions.size(); j++) {
            cout << "\t" << j << ":\t" << dat[i].labels[j] << endl;
            // cout << "\t" << j << ":\t" << dat[i].regions[j].first << " " << dat[i].regions[j].second << endl;
          }
        }
      }
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
      saveData(dat);
      return;
    }
  }

  bool isJoined = dat[0].regions.size()>1;

  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  if (!isJoined)
    for (auto &d : dat) {
      if (d.len > maxlen)
        maxlen = d.len;
      if (d.len < minlen)
        minlen = d.len;
    }
  else
    for (auto &it : dat[0].regions) {
      size_t len = it.second;
      if (len > maxlen)
        maxlen = len;
      if (len < minlen)
        minlen = len;
    }

  // adapt window size and interval
  size_t w = args.w;
  bool globalmode = w==0;
  size_t k = args.k;
  if (w == 0)
    w = minlen;       // default window = whole (smallest) seq.
  w = min(w, minlen); // biggest window = whole (smallest) seq.
  if (k == 0)
    k = max((size_t)1, w / 10); // default interval = 1/10 of window
  k = min(k, w);                // biggest interval = window size

  // array for results for all sequences in file
  size_t entries = numEntries(maxlen, w, k);
  size_t numMetrics = args.m == 'b' ? 2 : 1;
  vector<vector<double>> ys(numMetrics * max(dat.size(),dat[0].regions.size()), vector<double>(entries));

  int col = 0;
  if (!isJoined)
    for (auto &seq : dat) {
      size_t currw = globalmode ? seq.len : w;
      size_t currk = globalmode ? seq.len : k;
      if (args.m != 'r') {
        tick();
        mlComplexity(0, seq.len, currw, currk, ys[col++], seq);
        tock("mlComplexity");
      }

      if (args.m != 'm') {
        tick();
        runComplexity(0, seq.len, currw, currk, ys[col++], seq, true);
        tock("runComplexity");
      }
    }
  else
    for (auto &it : dat[0].regions) {
      size_t currw = globalmode ? it.second : w;
      size_t currk = globalmode ? it.second : k;
      if (args.m != 'r') {
        tick();
        mlComplexity(it.first, it.second, currw, currk, ys[col++], dat[0]);
        tock("mlComplexity");
      }

      if (args.m != 'm') {
        tick();
        runComplexity(it.first, it.second, currw, currk, ys[col++], dat[0], true);
        tock("runComplexity");
      }
    }

  if (!args.p) {
    cout << fixed << setprecision(4);
    if (args.g) { // print to be directly piped into some plot tool
      if (args.gf == 2) {
        for (size_t i = 0; i < ys.size(); i++) {
          for (size_t j = 0; j < entries; j++)
            cout << (j * k + w / 2) << "\t" << ys[i][j] << endl;
          cout << endl;
        }
      } else if (args.gf == 1) {
        cout << "DNALC_PLOT" << endl; //magic keyword
        gnuplotCode(w, k, ys.size()); // gnuplot control code
        // print column header (for plot labels)
        cout << "offset\t";
        for (size_t j = 0; j < max(dat.size(),dat[0].regions.size()); j++) { // columns for each seq
          string name = isJoined ? dat[0].labels[j] : dat[j].name;
          if (args.m != 'r')
            cout << "\"" << name << " (MC)\"\t";
          if (args.m != 'm')
            cout << "\"" << name << " (RC)\"\t";
        }
        cout << endl;
        printPlot(w, k, entries, ys); // print plot itself
      }
    } else { // just print resulting data
      printPlot(w, k, entries, ys);
    }
  }
}

int main(int argc, char *argv[]) {
  args.parse(argc, argv);

  // process files (or stdin, if none given)
  tick();
  if (!args.num_files) {
    processFile(nullptr);
  } else {
    for (size_t i = 0; i < args.num_files; i++) {
      processFile(args.files[i]);
    }
  }
  tock("total time");
}
