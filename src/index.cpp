#include <cassert>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
#include <algorithm>
using namespace std;

#include "args.h"  //args.p
#include "bench.h" //tick tock
#include "matchlength.h" //computeMLFact

#include "index.h"
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

  binwrite(o, (size_t)cd.regions.size());
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

  binwrite(o, (size_t)cd.fstRegionFact.size());
  for (auto i : cd.fstRegionFact)
    binwrite(o, i);
  binwrite(o, (size_t)cd.mlf.size());
  for (auto i : cd.mlf)
    binwrite(o, i);
}

const string magicstr = "BINIDX";

// serialize a series of sequences
bool saveData(vector<ComplexityData> &vec, char const *file) {
  bool joinedSeq = vec.size() == 1 && vec[0].regions.size() > 1;
  return with_file_out(file, [&](ostream &o) {
    for (auto c : magicstr) //magic sequence
      binwrite(o, c);
    binwrite(o, (char)joinedSeq);
    binwrite(o, (size_t)vec.size());
    for (auto &d : vec)
      serialize(o, d);
    return true;
  });
}

// load precomputed data from stdin (when file=nullptr) or some file
bool loadData(vector<ComplexityData> &cplx, char const *file, bool onlyInfo, size_t idx) {
  if (!file) {
    cerr << "ERROR: Can not load binary index file from pipe!"
      << " Please pass it as argument!" << endl;
    return false;
  }
  // return with_file(file, [&](istream &fin) {
  return with_mmap(file, [&](MMapReader &fin) {
    char tmp;
    for (size_t i = 0; i<magicstr.size(); i++) {
      binread(fin, tmp);
      if (tmp != magicstr[i]) {
        cerr << "ERROR: This does not look like an index file!" << endl;
        return false;
      }
    }

    char joinedSeq;
    binread(fin,joinedSeq);
    size_t n;
    binread(fin,n);
    for (size_t i = 0; i < n; i++) {
      //add new entry, THEN fill it (prevent copy)
      cplx.push_back(ComplexityData());
      auto &c = cplx.back();

      size_t namelen;
      binread(fin,namelen);
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

      size_t ffnum;
      binread(fin,ffnum);
      if (!onlyInfo)
        c.fstRegionFact.resize(ffnum);
      for (size_t j = 0; j < ffnum; j++) {
        size_t fi;
        binread(fin,fi);
        if (!onlyInfo)
          c.fstRegionFact[j] = fi;
      }

      size_t fnum;
      binread(fin,fnum);

      size_t start=0;
      size_t end=fnum;
      if (!onlyInfo) {
        if (!joinedSeq || !idx)
          c.mlf.resize(fnum);
        else {
          start = c.fstRegionFact[idx-1];
          end =  (size_t)idx==c.regions.size() ? fnum : c.fstRegionFact[idx];
          c.mlf.resize(end - start);
          fin.off += sizeof(size_t) * start;
        }
      }
      if (args.p)
        cerr << "Loading " << end-start << " of " << fnum << " factors..." << endl;
      for (size_t j = start; j < end; j++) {
        size_t fact;
        binread(fin,fact);
        if (!onlyInfo)
          c.mlf[j-start] = fact;
      }
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
    //first push, then fill (prevent copying in the end)
    cplx.push_back(ComplexityData());
    auto &c = cplx.back();

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

    size_t currreg=0;
    size_t idx=0;
    for (auto f : mlf.fact) {
      c.mlf.push_back(f);

      if (c.fstRegionFact.size() < c.regions.size() &&
          (c.fstRegionFact.empty() || (c.fstRegionFact.back()<f && f>=regions[currreg].first))) {
        c.fstRegionFact.push_back(idx);
        currreg++;
      }
      idx++;
    }

    if (args.p) {
      cout << seq.name << seq.comment << endl;
      // esa.print();
      cout << "ML-Factors for both strands (" << mlf.fact.size() << "):" << endl;
      mlf.print();
    }

    s.resize(s.size() / 2); // drop complementary seq.
    s.shrink_to_fit();
    esa.str = s.c_str();
    esa.n = s.size();
    /*
    tick();
    reduceEsa(esa);
    tock("reduceEsa");
    */

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
  }
}
