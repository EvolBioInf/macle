#include "util.h"
#include <cstring>
#include <cstdlib>

/* generate random sequence of given length */
char *randSeq(size_t n) {
  static char const *alphabet = "ACGT";
  char *s = new char[n + 2];
  s[n] = '$';
  s[n + 1] = '\0';
  for (size_t i = 0; i < n; i++)
    s[i] = alphabet[rand() % 4];
  return s;
}

/* fprintnf: print max of n characters of str onto fp;
 * add ... if str was truncated
 */
void fprintnf(FILE *fp, char const *str, int n) {
  int i, l, m;
  l = strlen(str);
  m = n < l ? n : l;
  for (i = 0; i < m; i++)
    fprintf(fp, "%c", str[i]);
  if (m < l)
    fprintf(fp, "...");
}
