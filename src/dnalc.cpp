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

typedef vector<pair<string,vector<double>>> ResultMat;

// given sequences from a fasta file, calculate match factors and runs
void extractData(vector<ComplexityData> &cplx, FastaFile &file, bool joinSeqs) {
  FastaFile *ff = &file;

  vector<string> labels;
  vector<pair<size_t,size_t>> regions;
  FastaFile virtfile;
  if (joinSeqs) { //construct concatenated sequence

    //extract region list from file if sequences are to be joined
    size_t offset=0;
    for (auto &it : file.seqs) {
      regions.push_back(make_pair(offset, it.seq.size()));
      labels.push_back(it.name+" "+it.comment);
      offset += it.seq.size();
    }

    string wholeseq = "";
    for (auto &seq : file.seqs) {
      wholeseq += seq.seq;
      seq.seq = ""; //free memory of separate sequences
    }
    virtfile.seqs.push_back(FastaSeq(file.filename, "", wholeseq));
    ff = &virtfile;
  }

  for (auto &seq : ff->seqs) {
    ComplexityData c;
    c.name = seq.name;
    c.gc = gcContent(seq.seq);
    c.len = seq.seq.size();
    c.labels = labels;
    c.regions = regions;

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

void indexInfo(vector<ComplexityData> const &dat) {
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
        cout << "\t" << (j+1) << ":\t" << dat[i].labels[j] << endl;
        // cout << dat[i].regions[j].first << " " << dat[i].regions[j].second << endl;
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

// This complicated function calculates the complexity data depending on mode.
// The input data can be "joined" -> one single sequence with "regions", or all
// sequences in the input are separate.
// Also the user can use global mode (one value per complexity and sequence),
// and the user can choose a region or sequence to process.
// If NOT joined: no chosen seqnum -> compute for all separately, otherwise only given sequence
// If joined: no chosen seqnum -> compute for complete sequence, otherwise only one region
ResultMat calcComplexities(size_t &w, size_t &k, char m, size_t seqnum, vector<ComplexityData> const &dat) {
  bool globalMode = w==0;   // output one number (window = whole sequence)?
  bool isJoined = dat[0].regions.size()>1; //which kind is the input data?
  size_t containedSeqs = max(dat[0].labels.size(), dat.size()); //how many (sub-)sequences are there?

  size_t maxlen = 0;        // max implies the domain of the plot
  size_t minlen = SIZE_MAX; // min restricts the reasonable window sizes
  if (seqnum || isJoined)
    maxlen = minlen = seqnum && isJoined ? dat[0].regions[seqnum-1].second : dat[0].len;
  else //we need to consider the separate sequences, non-joined mode
    for (auto &d : dat) {
      maxlen = max(maxlen, d.len);
      minlen = min(minlen, d.len);
    }

  // adapt window size and interval
  if (w == 0)
    w = minlen;       // default window = whole (smallest) seq.
  w = min(w, minlen); // biggest window = whole (smallest) seq.
  if (k == 0)
    k = max((size_t)1, w / 10); // default interval = 1/10 of window
  k = min(k, w);                // biggest interval = window size

  // array for results for all sequences in file
  size_t entries = numEntries(maxlen, w, k);
  size_t numMetrics = m == 'b' ? 2 : 1; //complexity arrays per seq.
  size_t usedSeqs = seqnum || isJoined ? 1 : containedSeqs;
  ResultMat ys(numMetrics * usedSeqs, make_pair("", vector<double>(entries)));

  int col = 0;
  int start = seqnum ? max(seqnum-1,(size_t)0)       : 0;
  int end   = seqnum ? min(seqnum-1,containedSeqs-1) : (isJoined ? start : containedSeqs-1);
  for (int i=start; i<=end; i++) {
    auto &data    = isJoined ? dat[0]                                           : dat[i];
    string name   = isJoined ? (seqnum ? dat[0].labels[i] : dat[0].name)        : dat[i].name;
    size_t offset = isJoined ? (seqnum ? dat[0].regions[i].first : 0)           : 0;
    size_t len    = isJoined ? (seqnum ? dat[0].regions[i].second : dat[0].len) : dat[i].len;
    size_t currw  = globalMode ? len : w;
    size_t currk  = globalMode ? len : k;

    if (m != 'r') {
      tick();
      ys[col].first = name + " (MC)";
      mlComplexity(offset, len, currw, currk, ys[col].second, data);
      tock("mlComplexity");
      col++;
    }
    if (m != 'm') {
      tick();
      ys[col].first = name + " (RC)";
      runComplexity(offset, len, currw, currk, ys[col].second, data, true);
      tock("runComplexity");
      col++;
    }
  }
  return ys;
}

//load / extract data, show results
void processFile(char const *file) {
  vector<ComplexityData> dat;
  if (args.i) { //load from index
    tick();
    if (!loadData(dat, file))
      return;
    tock("loadData");
    if (args.l) { //list index file contents and exit
      indexInfo(dat);
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
