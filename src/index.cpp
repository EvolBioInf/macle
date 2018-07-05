#include <cassert>
#include <iostream>
#include <fstream>
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

const string magicstr = "BINIDX";

bool saveData(ComplexityData &cd, char const *file) {
  assert(cd.regions.size() == cd.labels.size());
  return with_file_out(file, [&](ostream &o) {
    for (auto c : magicstr) //magic sequence
      binwrite(o, c);

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
    for (auto &l : cd.labels) {
      binwrite(o, l.size());
      for (size_t i=0; i<MAX_LABEL_LEN; i++)
        binwrite(o, i < l.size() ? l[i] : (char)0);
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

    return true;
  });
}

//take stream, consume magic string, return success value
bool readMagic(istream &fin) {
  char tmp;
  for (size_t i = 0; i<magicstr.size(); i++) {
    binread(fin, tmp);
    if (tmp != magicstr[i]) {
      return false;
    }
  }
  return true;
}

// load precomputed data from stdin (when file=nullptr) or some file
bool loadData(ComplexityData &dat, char const *file, bool onlyInfo) {
  if (!file) {
    cerr << "ERROR: Can not load binary index file from pipe!"
      << " Please pass it as argument!" << endl;
    return false;
  }
  return with_file_in(file, [&](istream &fin) {
  // return with_mmap(file, [&](MMapReader &fin) {
    if (!readMagic(fin)) {
      cerr << "ERROR: This does not look like an index file!" << endl;
      return false;
    }

    char tmp;
    size_t namelen;
    binread(fin,namelen);
    for (size_t j=0; j<namelen; j++) {
      binread(fin,tmp);
      dat.name += tmp;
    }
    binread(fin,dat.len);
    binread(fin,dat.gc);

    size_t rnum;
    binread(fin,rnum);
    dat.regions.resize(rnum);
    dat.labels.resize(rnum);
    for (size_t j = 0; j < rnum; j++) {
      size_t s, l;
      binread(fin,s);
      binread(fin,l);
      dat.regions[j] = make_pair(s, l);
    }
    for (size_t j = 0; j < rnum; j++) {
      size_t lbllen;
      binread(fin,lbllen);
      for (size_t k=0; k<MAX_LABEL_LEN; k++) {
        binread(fin,tmp);
        if (k < lbllen)
          dat.labels[j] += tmp;
      }
    }

    binread(fin, dat.numbad);

    size_t bnum;
    binread(fin,bnum);
    if (!onlyInfo)
      dat.bad.resize(bnum);
    for (size_t j = 0; j < bnum; j++) {
      size_t l, r;
      binread(fin,l);
      binread(fin,r);
      if (!onlyInfo)
        dat.bad[j] = make_pair(l, r);
    }

    size_t ffnum;
    binread(fin,ffnum);
    if (!onlyInfo)
      dat.fstRegionFact.resize(ffnum);
    for (size_t j = 0; j < ffnum; j++) {
      size_t fi;
      binread(fin,fi);
      if (!onlyInfo)
        dat.fstRegionFact[j] = fi;
    }

    size_t fnum;
    binread(fin,fnum);

    if (!onlyInfo)
      dat.mlf.resize(fnum);
    for (size_t j = 0; j < fnum; j++) {
      size_t fact;
      binread(fin,fact);
      if (!onlyInfo)
        dat.mlf[j] = fact;
    }

    return true;
  }, ios::in|ios::binary);
}

bool renameRegions(char const *file, vector<string> const &names) {
  if (!file) {
    cerr << "ERROR: Can not load binary index file from pipe!"
      << " Please pass it as argument!" << endl;
    return false;
  }
  return with_file(file, [&](fstream &fs) {
    char tmp;
    for (size_t i = 0; i<magicstr.size(); i++) {
      binread(fs, tmp);
      if (tmp != magicstr[i]) {
        cerr << "ERROR: This does not look like an index file!" << endl;
        return false;
      }
    }
    size_t tmpsz;
    binread(fs,tmpsz);
    for (size_t j=0; j<tmpsz; j++)
      binread(fs,tmp);
    double gc;
    binread(fs,tmpsz);
    binread(fs,gc);

    size_t rnum;
    binread(fs,rnum);
    if (rnum != names.size()) {
      cerr << "ERROR: number of given names and regions does not match!" << endl;
      return false;
    }
    for (size_t j = 0; j < rnum; j++) {
      binread(fs,tmpsz);
      binread(fs,tmpsz);
    }
    fs.seekp(fs.tellg());
    for (auto &s : names) {
      binwrite(fs, s.size());
      for (size_t i=0; i<MAX_LABEL_LEN; i++)
        binwrite(fs, i < s.size() ? s[i] : (char)0);
    }
    return true;
  }, fstream::in|fstream::out|fstream::binary);
}

// given sequences from a fasta file, calculate match factors and runs
void extractData(ComplexityData &dat, FastaFile &file) {
  //construct concatenated sequence:
  string s = "";
  size_t offset=0;
  for (auto &it : file.seqs) { //extract region list from file
    dat.regions.push_back(make_pair(offset, it.seq.size()));
    dat.labels.push_back(it.name /* +" "+it.comment */);
    offset += it.seq.size();

    s += it.seq;
    it.seq = ""; //free memory of separate sequences
  }

  dat.name = file.filename;
  dat.gc = gcContent(s);
  dat.len = s.size();

  s = s + "$" + revComp(s) + "$";
  tick();
  Esa esa(s.c_str(), s.size()); // esa for seq+$+revseq+$
  tock("getEsa (both strands)");

  tick();
  Fact mlf;
  computeMLFact(mlf, esa);
  tock("computeMLFact");

  s.resize(s.size() / 2); // drop complementary seq.
  s.shrink_to_fit();
  mlf.str = esa.str = s.c_str();
  esa.n = s.size();

  size_t currreg=0;
  size_t idx=0;
  for (auto f : mlf.fact) {
    dat.mlf.push_back(f);
    if (dat.fstRegionFact.size() < dat.regions.size() &&
        (dat.fstRegionFact.empty() || (dat.fstRegionFact.back()<f && f>=dat.regions[currreg].first))) {
      dat.fstRegionFact.push_back(idx);
      currreg++;
    }
    idx++;
  }

  if (args.p) {
    // esa.print();
    cout << "ML-Factors on first strand (" << mlf.fact.size() << "):" << endl;
    mlf.print();
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
      dat.bad.push_back(make_pair(start, i - 1));
    } else if (!insidebad && !valid) {
      start = i;
      insidebad = true;
    }
  }
  if (insidebad) // push last one, if we are inside
    dat.bad.push_back(make_pair(start, s.size() - 1));

  //calculate total number of bad nucleotides for global mode
  dat.numbad=0;
  for (auto &bad : dat.bad)
    dat.numbad += bad.second - bad.first + 1;

  tock("find bad intervals");
}
