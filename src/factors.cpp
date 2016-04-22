#include "prelude.h"
#include "factors.h"

void freeFact(Fact *f) {
  free(f->fact);
  if (f->lpf)
    free(f->lpf);
  if (f->prevOcc)
    free(f->prevOcc);
  free(f);
}

void printFact(Fact *mlf) {
  char tmp;
  char* str = (char*)mlf->str; //we write there, but restore it back!
  for (size_t i = 0; i < mlf->n; i++) {
    size_t start = mlf->fact[i];
    size_t end = i < mlf->n - 1 ? mlf->fact[i + 1] : mlf->strLen;
    tmp = mlf->str[end];
    str[end] = '\0';
    printf("%s%s", &(mlf->str[start]), i < mlf->n - 1 ? "." : "\n");
    str[end] = tmp;
  }
}

extern inline size_t factLen(Fact *f, size_t i);
