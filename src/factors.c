#include "prelude.h"
#include "factors.h"

void freeFact(Fact *f) {
  free(f->fact);
  if (f->lpf)
    free(f->lpf);
  free(f);
}

void printFact(Fact *mlf) {
  char tmp;
  for (size_t i = 0; i < mlf->n; i++) {
    size_t start = mlf->fact[i];
    size_t end = i < mlf->n - 1 ? mlf->fact[i + 1] : mlf->strLen;
    tmp = mlf->str[end];
    mlf->str[end] = '\0';
    printf("%s%s", &(mlf->str[start]), i < mlf->n - 1 ? "." : "\n");
    mlf->str[end] = tmp;
  }
}
