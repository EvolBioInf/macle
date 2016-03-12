#include "minunit.h"
#include <string.h>
#include <time.h>

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
  Esa *esa = getEsa(s, n); // calculate esa, including $
  Fact *lzf = computeLZFact(esa, false);

  size_t plen;
  Periodicity *ps = getPeriodicities2(esa, &plen);
  mu_assert_eq(8, plen, "wrong number of periodicities detected (easy algorithm)");
  free(ps);

  List **pl = getPeriodicityLists(false, lzf, esa, &plen);
  mu_assert_eq(8, plen, "wrong number of periodicities detected");
  freePeriodicityLists(pl, n);

  freeFact(lzf);
  freeEsa(esa);
  return NULL;
}

char *test_onlyRuns() {
  char *s = "AAAAAAAAGCGCGCGCGCGCGCGTTTTTTTTTTTTACTACTACTACTACTACTA$";
  size_t n = strlen(s);
  Esa *esa = getEsa(s, n); // calculate esa, including $
  Fact *lzf = computeLZFact(esa, false);

  size_t plen;
  Periodicity *ps = getPeriodicities2(esa, &plen);
  mu_assert_eq(16, plen, "wrong number of periodicities detected (easy algorithm)");
  free(ps);

  ps = getPeriodicities(false, lzf, esa, &plen);
  mu_assert_eq(16, plen, "wrong number of periodicities detected");
  free(ps);

  ps = getPeriodicities(true, lzf, esa, &plen);
  mu_assert_eq(4, plen, "wrong number of runs detected");
  free(ps);

  freeFact(lzf);
  freeEsa(esa);
  return NULL;
}

char *test_randomSequence() {
  size_t n = 1000000;
  char *s = randSeq(n);
  fprintnf(stdout, s, 80);
  printf("\n");
  Esa *esa = getEsa(s, n + 1); // calculate esa, including $
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

      fprintnf(stdout, s + ps2[i + missed].b - 1, 80);
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
  free(s);
  mu_assert_eq(plen2, plen, "number of periodicities does not match");
  return NULL;
}

// compare: getting all, then filtering afterwards (as proposed in book)
// vs. getting only runs directly (for better performance)
char *test_randomOnlyRuns() {
  size_t n = 1000000;
  char *s = randSeq(n);
  n++; //$ border
  fprintnf(stdout, s, 80);
  printf("\n");
  Esa *esa = getEsa(s, n); // calculate esa
  Fact *lzf = computeLZFact(esa, false);

  size_t plen, plen2;
  // get all, then remove non-runs
  List **pl = getPeriodicityLists(false, lzf, esa, &plen);
  for (size_t i = 0; i < lzf->strLen; i++) {
    Periodicity *lastp = NULL;
    List *last = NULL;
    for (eachListItem(curr, pl[i])) {
      Periodicity *currp = (Periodicity *)curr->value;
      if (lastp && lastp->b == currp->b && lastp->e == currp->e) {
        last->next = curr->next;
        free(currp);
        free(curr);
        plen--;
        curr = last;
      } else {
        last = curr;
        lastp = currp;
      }
    }
  }
  freePeriodicityLists(pl, n);

  // get just runs
  Periodicity *ps = getPeriodicities(true, lzf, esa, &plen2);
  free(ps);

  freeFact(lzf);
  freeEsa(esa);
  free(s);
  mu_assert_eq(plen2, plen, "number of runs does not match");
  return NULL;
}

char *all_tests() {
  srand(time(NULL));
  mu_suite_start();
  mu_run_test(test_knownExample);
  mu_run_test(test_onlyRuns);
  for (size_t i = 0; i < 3; i++)
    mu_run_test(test_randomSequence);
  for (size_t i = 0; i < 3; i++)
    mu_run_test(test_randomOnlyRuns);
  return NULL;
}
RUN_TESTS(all_tests)
