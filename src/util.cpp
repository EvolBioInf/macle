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
  uniform_int_distribution<int> distr(0, 3);
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
