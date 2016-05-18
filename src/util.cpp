#include <cstdio>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstring>

//for open_or_fail
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

#include "lempelziv.h"
#include "util.h"
#include "args.h"
using namespace std;

// generate random DNA seq of given length
string randSeq(size_t n, double gc) {
  string s;
  s.resize(n);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine gen(seed);
  bernoulli_distribution coin(0.5);
  uniform_real_distribution<double> rdistr(0, 1.0);
  for (size_t i = 0; i < n; i++) {
    if (rdistr(gen)<gc)
      s[i] = coin(gen) ? 'G' : 'C';
    else
      s[i] = coin(gen) ? 'A' : 'T';
  }
  return s;
}

// generate random sequence of given length
string randSeq(size_t n, string alphabet) {
  string s;
  s.resize(n);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine gen(seed);
  uniform_int_distribution<int> distr(0, alphabet.size() - 1);
  for (size_t i = 0; i < n; i++)
    s[i] = alphabet[distr(gen)];
  return s;
}

string randRun(size_t n, size_t l, string alphabet) {
  if (l > n)
    l = n;
  string s;
  string rep = randSeq(l, alphabet);
  while (s.size() < n)
    s += rep;
  s.resize(n);
  return s;
}

// calculate the GC content
double gcContent(string const &s) {
  size_t gc = 0;
  for (auto c : s)
    if (c == 'g' || c == 'G' || c == 'c' || c == 'C')
      gc++;
  return (double)gc / s.size();
}

// returns reverse complement DNA string
string revComp(string const &s) {
  string r = s;
  for (char &c : r)
    switch (c) {
    case 'A':
      c = 'T';
      break;
    case 'T':
      c = 'A';
      break;
    case 'G':
      c = 'C';
      break;
    case 'C':
      c = 'G';
      break;
    }
  reverse(r.begin(), r.end());
  return r;
}

/* open_or_fail: open file on system level and report on error */
int open_or_fail(char const *fname, int flag) {
  int fd = open(fname, flag, 0);
  if (fd < 0) {
    fprintf(stderr, "open(%s, %i) failed: %s\n", fname, flag, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return fd;
}

/* fopen_or_fail: open file on system level and report on error */
FILE *fopen_or_fail(char const *fname, char const *flags) {
  FILE *f = fopen(fname, flags);
  if (!f) {
    fprintf(stderr, "fopen(%s, %s) failed: %s\n", fname, flags, strerror(errno));
    exit(EXIT_FAILURE);
  }
  return f;
}

// fprintnf: print max of n characters of str onto fp;
// add ... if str was truncated
void fprintnf(FILE *fp, char const *str, int n) {
  int i, l, m;
  l = strlen(str);
  m = n < l ? n : l;
  for (i = 0; i < m; i++)
    fprintf(fp, "%c", str[i]);
  if (m < l)
    fprintf(fp, "...");
}

vector<list<Periodicity>> getRuns(string const &seq) {
  bool b = args.b;
  args.b = false;
  Esa esa(seq.c_str(), seq.size()); // esa for sequence+$
  Fact lzf;
  computeLZFact(lzf, esa, false);
  size_t pnum;
  auto ls=getPeriodicityLists(true, lzf, esa, pnum);
  args.b = b;
  return ls;
}
