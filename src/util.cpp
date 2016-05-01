#include <cstdlib>
#include <string>
#include <cstring>
using namespace std;

// generate random sequence of given length
string randSeq(size_t n) {
  static char const *alphabet = "ACGT";
  string s;
  s.reserve(n + 2);
  for (size_t i = 0; i < n; i++)
    s += alphabet[rand() % 4];
  s += "$";
  return s;
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
