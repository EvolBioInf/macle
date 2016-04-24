#include <cstdio>
#include "factors.h"

Fact::~Fact() {
  delete[] this->fact;
  if (this->lpf)
    delete[] this->lpf;
  if (this->prevOcc)
    delete[] this->prevOcc;
}

void Fact::print() const {
  char tmp;
  char *s = (char *)this->str; // we write there, but restore it back!
  for (size_t i = 0; i < this->n; i++) {
    size_t start = this->fact[i];
    size_t end = i < this->n - 1 ? this->fact[i + 1] : this->strLen;
    tmp = this->str[end];
    s[end] = '\0';
    printf("%s%s", &(this->str[start]), i < this->n - 1 ? "." : "\n");
    s[end] = tmp;
  }
}

extern inline size_t factLen(Fact &f, size_t i);
