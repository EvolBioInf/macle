#include <cassert>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
using namespace std;

#include "args.h"  //args.p
#include "bench.h" //tick tock
#include "matchlength.h" //computeMLFact
#include "lempelziv.h"   //computeLZFact

#include "index.h"
#include "periodicity.h"
#include "util.h"

template<typename T> void binwrite(ostream &o, T x) {
  o.write(reinterpret_cast<char*>(&x),sizeof(x));
}
template<typename T> void binread(istream &i, T &x) {
  i.read(reinterpret_cast<char*>(&x),sizeof(x));
}
template<typename T> void binget(MMapReader &i, size_t offbytes, T &x) {
  x = *reinterpret_cast<T*>(i.dat+i.off+offbytes);
}
template<typename T> void binread(MMapReader &i, T &x) {
  binget(i, 0, x);
  i.off += sizeof(x);
}

// serialize relevant data of a sequence
void serialize(ostream &o, ComplexityData const &cd) {
  assert(cd.regions.size() == cd.labels.size());

  binwrite(o, (size_t)cd.name.size());
  for (auto c : cd.name)
    binwrite(o, c);
  binwrite(o, cd.len);
  binwrite(o, cd.gc);

  binwrite(o, cd.regions.size());
  for (auto i : cd.regions){
    binwrite(o, i.first);
    binwrite(o, i.second);
  }
  for (auto &l : cd.labels){
    binwrite(o, l.size());
    for (char c : l)
      binwrite(o, c);
  }

  binwrite(o, cd.numbad);
  binwrite(o, (size_t)cd.bad.size());
  for (auto i : cd.bad){
    binwrite(o, i.first);
    binwrite(o, i.second);
  }

  binwrite(o, (size_t)cd.mlf.size());
  for (auto i : cd.mlf)
    binwrite(o, i);
  size_t pernum = 0;
  for (auto &l : cd.pl)
    pernum += l.size();
  binwrite(o, pernum);
  for (auto &l : cd.pl)
    for (auto p : l) {
      binwrite(o, p.b);
      binwrite(o, p.e);
      binwrite(o, p.l);
    }
}

// serialize a series of sequences
bool saveData(vector<ComplexityData> &vec, char const *file) {
  return with_file_out(file, [&](ostream &o) {
    binwrite(o, (size_t)vec.size());
    for (auto &d : vec)
      serialize(o, d);
    return true;
  });
}

// load precomputed data from stdin (when file=nullptr) or some file
bool loadData(vector<ComplexityData> &cplx, char const *file, bool onlyInfo) {
  if (!file) {
    cerr << "ERROR: Can not load binary index file from pipe!"
      << " Please pass it as argument!" << endl;
    return false;
  }
  // return with_file(file, [&](istream &fin) {
  return with_mmap(file, [&](MMapReader &fin) {
    size_t n;
    binread(fin,n);
    for (size_t i = 0; i < n; i++) {
      ComplexityData c;
      size_t namelen;
      binread(fin,namelen);
      char tmp;
      for (size_t j=0; j<namelen; j++) {
        binread(fin,tmp);
        c.name += tmp;
      }
      binread(fin,c.len);
      binread(fin,c.gc);

      size_t rnum;
      binread(fin,rnum);
      c.regions.resize(rnum);
      c.labels.resize(rnum);
      for (size_t j = 0; j < rnum; j++) {
        size_t s, l;
        binread(fin,s);
        binread(fin,l);
        c.regions[j] = make_pair(s, l);
      }
      for (size_t j = 0; j < rnum; j++) {
        size_t lbllen;
        binread(fin,lbllen);
        for (size_t k=0; k<lbllen; k++) {
          binread(fin,tmp);
          c.labels[j] += tmp;
        }
      }

      binread(fin, c.numbad);

      size_t bnum;
      binread(fin,bnum);
      if (!onlyInfo)
        c.bad.resize(bnum);
      for (size_t j = 0; j < bnum; j++) {
        size_t l, r;
        binread(fin,l);
        binread(fin,r);
        if (!onlyInfo)
          c.bad[j] = make_pair(l, r);
      }

      size_t fnum;
      binread(fin,fnum);
      if (!onlyInfo)
        c.mlf.resize(fnum);
      for (size_t j = 0; j < fnum; j++) {
        size_t fact;
        binread(fin,fact);
        if (!onlyInfo)
          c.mlf[j] = fact;
      }

      size_t pnum;
      binread(fin,pnum);
      if (!onlyInfo)
        c.pl.resize(c.len + 2); // +2 because of $ and one more for safety
      for (size_t j = 0; j < pnum; j++) {
        size_t b, e, l;
        binread(fin,b);
        binread(fin,e);
        binread(fin,l);
        if (!onlyInfo)
          c.pl[b].push_back(Periodicity(b, e, l));
      }

      cplx.push_back(c);
    }
    return true;
  }); //, ios::binary);
}

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
