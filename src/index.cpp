#include <cassert>
#include <iostream>
#include <vector>
#include <utility>
#include <functional>
using namespace std;

#include "index.h"
#include "periodicity.h"
#include "util.h"
using namespace std;

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

  binwrite(o, (size_t)cd.bad.size());
  for (auto i : cd.bad){
    binwrite(o, i.first);
    binwrite(o, i.second);
  }
  binwrite(o, cd.numbad);
}

// serialize a series of sequences
bool saveData(vector<ComplexityData> &vec, char const *file) {
  return with_file(file, [&](ostream &o) {
    binwrite(o, (size_t)vec.size());
    for (auto &d : vec)
      serialize(o, d);
    return true;
  });
}

// load precomputed data from stdin (when file=nullptr) or some file
bool loadData(vector<ComplexityData> &cplx, char const *file) {
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

      size_t fnum;
      binread(fin,fnum);
      c.mlf.resize(fnum);
      for (size_t j = 0; j < fnum; j++) {
        binread(fin,c.mlf[j]);
      }

      size_t pnum;
      binread(fin,pnum);
      c.pl.resize(c.len);
      for (size_t j = 0; j < pnum; j++) {
        size_t b, e, l;
        binread(fin,b);
        binread(fin,e);
        binread(fin,l);
        c.pl[b].push_back(Periodicity(b, e, l));
      }

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

      size_t bnum;
      binread(fin,bnum);
      c.bad.resize(bnum);
      for (size_t j = 0; j < bnum; j++) {
        size_t l, r;
        binread(fin,l);
        binread(fin,r);
        c.bad[j] = make_pair(l, r);
      }
      binread(fin, c.numbad);

      cplx.push_back(c);
    }
    return true;
  }); //, ios::binary);
}

