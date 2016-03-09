#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"
#include "lempelziv.h"
#include "periodicity.h"
#include "stringUtil.h"

// helper: compare periodicities
int cmpPer(const void *a, const void *b) {
  Periodicity *p1 = (Periodicity *)a;
  Periodicity *p2 = (Periodicity *)b;
  int db = p1->b - p2->b;
  int de = p1->e - p2->e;
  int dl = p1->l - p2->l;
  if (db)
    return db;
  if (de)
    return de;
  return dl;
}

char *test_knownExample() {
  char *s = "AACCAACCAACCAA$"; // from Ohlebusch book
  size_t n = strlen(s);
  Esa *esa = getEsa(s, n + 1); // calculate esa, including $
  Fact *lzf = computeLZFact(esa, false);

  size_t plen;
  Periodicity *ps = getPeriodicities(false, lzf, esa, &plen);
  mu_assert(plen == 8, "wrong number of periodicities detected");

  free(ps);
  freeFact(lzf);
  freeEsa(esa);
  return NULL;
}

char *test_onlyRuns() {
  char *s = "AAAAAAAAGCGCGCGCGCGCGCGTTTTTTTTTTTTACTACTACTACTACTACTA$";
  size_t n = strlen(s);
  Esa *esa = getEsa(s, n + 1); // calculate esa, including $
  Fact *lzf = computeLZFact(esa, false);

  size_t plen;
  Periodicity *ps = getPeriodicities(false, lzf, esa, &plen);
  mu_assert(plen == 16, "wrong number of periodicities detected");
  free(ps);

  ps = getPeriodicities2(esa, &plen);
  mu_assert(plen == 16, "wrong number of periodicities detected (easy algorithm)");
  free(ps);

  ps = getPeriodicities(true, lzf, esa, &plen);
  mu_assert(plen == 4, "wrong number of runs detected");
  free(ps);

  freeFact(lzf);
  freeEsa(esa);
  return NULL;
}

char *test_compareWithReference() {
  Sequence *seq = readFastaFromFile("Data/broken.fa");
  Esa *esa = getEsa(seqStr(seq, 0), seqLen(seq, 0) + 1); // calculate esa, including $
  Fact *lzf = computeLZFact(esa, false);

  size_t plen, plen2;
  Periodicity *ps, *ps2;

  ps = getPeriodicities(false, lzf, esa, &plen);
  qsort(ps, plen, sizeof(Periodicity), cmpPer);
  ps2 = getPeriodicities2(esa, &plen2);
  qsort(ps2, plen2, sizeof(Periodicity), cmpPer);

  size_t missed = 0;
  for (size_t i = 0; i + missed < plen; i++) {
    if (ps[i].b != ps2[i + missed].b || ps[i].e != ps2[i + missed].e ||
        ps[i].l != ps2[i + missed].l) {
      printf("found: ");
      printPeriodicity(&ps[i]);
      printf("should: ");
      printPeriodicity(&ps2[i + missed]);
      fprintnf(stdout, seq->seq + ps2[i + missed].b - 1, 80);
      printf("\nrelevant LZ-Factor:\n");
      size_t fact = 0;
      for (size_t j = 0; j < lzf->n; j++) {
        size_t b = lzf->fact[j];
        size_t e = lzf->fact[j] + factLen(lzf, j);
        if (b <= ps2[i + missed].b - 1 && e >= ps2[i + missed].e - 1) {
          printf("%zu: ", j);
          fprintnf(stdout, lzf->str + lzf->fact[j], factLen(lzf, j));
          printf("\n");
          fact = lzf->fact[j];
        }
      }
      int64_t prev = fact;
      fact = 0;
      while (prev != -1) {
        fact = prev;
        prev = lzf->prevOcc[fact];
        printf("pos: %ld\n", prev);
      }

      missed++;
    }
  }

  free(ps);
  free(ps2);
  freeFact(lzf);
  freeEsa(esa);
  freeSequence(seq);
  mu_assert(plen == plen2, "periodicities do not match");
  return NULL;
}

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_knownExample);
  mu_run_test(test_onlyRuns);
  mu_run_test(test_compareWithReference);
  return NULL;
}
RUN_TESTS(all_tests)
