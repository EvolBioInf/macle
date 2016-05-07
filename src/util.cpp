#include <cstdlib>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <cstring>
#include "util.h"
using namespace std;

// generate random sequence of given length
string randSeq(size_t n, string alphabet) {
  string s;
  s.resize(n);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine gen(seed);
  uniform_int_distribution<int> distr(0,3);
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
  for (auto it = s.begin(); it != s.end(); it++)
    if (*it == 'g' || *it == 'G' || *it == 'c' || *it == 'C')
      gc++;
  return (double)gc / s.size();
}

// returns reverse complement DNA string
string revComp(string const &s) {
  string r = s;
  for (auto it = r.begin(); it != r.end(); it++)
    switch (*it) {
    case 'A':
      *it = 'T';
      break;
    case 'T':
      *it = 'A';
      break;
    case 'G':
      *it = 'C';
      break;
    case 'C':
      *it = 'G';
      break;
    }
  reverse(r.begin(), r.end());
  return r;
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
