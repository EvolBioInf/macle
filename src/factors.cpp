#include <iostream>
using namespace std;

#include "factors.h"

void Fact::print() const {
  char tmp;
  size_t n = fact.size();
  char *s = (char *)this->str; // we write there, but restore it back!
  for (size_t i = 0; i < n; i++) {
    size_t start = this->fact[i];
    size_t end = i < n - 1 ? this->fact[i + 1] : this->strLen;
    tmp = this->str[end];
    s[end] = '\0';
    cout << this->str+start << (i < n - 1 ? "." : "\n");
    s[end] = tmp;
  }
}

extern inline size_t factLen(Fact &f, size_t i);
