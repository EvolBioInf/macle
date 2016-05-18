#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <functional>
using namespace std;

#include <cstdio>

#include "index.h"
#include "periodicity.h"
#include "util.h"
using namespace std;

template<typename T> void binwrite(T x) {
  cout.write(reinterpret_cast<const char*>(&x),sizeof(x));
}
template<typename T> void binread(istream &i, T &x) {
  i.read(reinterpret_cast<char*>(&x),sizeof(x));
}
bool with_file(char const *file, function<bool(istream&)> lambda, ios_base::openmode mode=ios_base::in) {
  istream *finP = &cin;
  ifstream fs;
  if (file) {
    // fs = ifstream(file); //does not work with older g++?
    fs.open(file, mode);
    if (!fs.is_open()) {
      cerr << "ERROR: Could not open file: " << file << endl;
      return false;
    }
    finP = &fs;
  }
  istream &fin = *finP;
  bool ret = lambda(fin);
  if (fs.is_open())
    fs.close();
  return ret;
}

// serialize relevant data of a sequence
void serialize(ComplexityData const &cd) {
  binwrite((size_t)cd.name.size());
  for (auto c : cd.name)
    binwrite(c);
  binwrite(cd.len);
  binwrite(cd.gc);
  binwrite((size_t)cd.mlf.size());
  for (auto i : cd.mlf)
    binwrite(i);
  size_t pernum = 0;
  for (auto &l : cd.pl)
    pernum += l.size();
  binwrite(pernum);
  for (auto &l : cd.pl)
    for (auto p : l) {
      binwrite(p.b);
      binwrite(p.e);
      binwrite(p.l);
    }
  binwrite((size_t)cd.bad.size());
  for (auto i : cd.bad){
    binwrite(i.first);
    binwrite(i.second);
  }
}

// serialize a series of sequences
void saveData(vector<ComplexityData> &vec) {
  binwrite((size_t)vec.size());
  for (auto &d : vec)
    serialize(d);
}

// load precomputed data from stdin (when file=nullptr) or some file
bool loadData(vector<ComplexityData> &cplx, char const *file) {
  return with_file(file, [&](istream &fin) {
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
      for (size_t j = 0; j < fnum; j++)
        binread(fin,c.mlf[j]);

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
      size_t bnum;
      binread(fin,bnum);
      for (size_t j = 0; j < bnum; j++) {
        size_t l, r;
        binread(fin,l);
        binread(fin,r);
        c.bad.push_back(make_pair(l, r));
      }

      cplx.push_back(c);
    }
    return true;
  }, ios::binary);
}

