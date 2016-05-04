#include "minunit.h"
#include "util.h"
#include <cstring>
#include <ctime>

#include "esa.h"
#include "lempelziv.h"
#include "periodicity.h"

#include <vector>
#include <list>
#include <algorithm>
using namespace std;

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

char const *test_knownExample() {
  char const *s = "AACCAACCAACCAA$"; // from Ohlebusch book
  size_t n = strlen(s);
  Esa esa(s, n); // calculate esa, including $
  Fact lzf;
  computeLZFact(lzf, esa, false);

  auto ps = getPeriodicities2(esa);
  mu_assert_eq(8, ps.size(), "wrong number of periodicities detected (easy algorithm)");

  size_t plen;
  auto pl = getPeriodicityLists(false, lzf, esa, plen);
  mu_assert_eq(8, plen, "wrong number of periodicities detected");

  return NULL;
}

char const *test_onlyRuns() {
  char const *s = "AAAAAAAAGCGCGCGCGCGCGCGTTTTTTTTTTTTACTACTACTACTACTACTA$";
  size_t n = strlen(s);
  Esa esa(s, n); // calculate esa, including $
  Fact lzf;
  computeLZFact(lzf, esa, false);

  auto ps2 = getPeriodicities2(esa);
  mu_assert_eq(16, ps2.size(), "wrong number of periodicities detected (easy algorithm)");

  size_t plen;
  auto ps = getPeriodicities(false, lzf, esa, plen);
  mu_assert_eq(16, plen, "wrong number of periodicities detected");

  ps = getPeriodicities(true, lzf, esa, plen);
  mu_assert_eq(4, plen, "wrong number of runs detected");

  return NULL;
}

char const *test_randomSequence() {
  size_t n = 1000000;
  string str = randSeq(n);
  char const *s = str.c_str();
  fprintnf(stdout, s, 80);
  printf("\n");
  Esa esa(s, n + 1); // calculate esa, including $
  Fact lzf;
  computeLZFact(lzf, esa, false);

  size_t plen;
  auto ps = getPeriodicities(false, lzf, esa, plen);
  sort(ps.begin(), ps.end(),
       [](Periodicity a, Periodicity b) { return cmpPer(&a, &b) < 0; });

  auto ps2 = getPeriodicities2(esa);
  sort(ps2.begin(), ps2.end(),
       [](Periodicity a, Periodicity b) { return cmpPer(&a, &b) < 0; });

  size_t missed = 0;
  for (size_t i = 0; i + missed < plen; i++) {
    if (ps[i].b != ps2[i + missed].b || ps[i].e != ps2[i + missed].e ||
        ps[i].l != ps2[i + missed].l) {
      printf("found: ");
      printPeriodicity(ps[i]);
      printf("should: ");
      printPeriodicity(ps2[i + missed]);

      fprintnf(stdout, s + ps2[i + missed].b - 1, 80);
      printf("\nrelevant LZ-Factor:\n");
      size_t fact = 0;
      for (size_t j = 0; j < lzf.fact.size(); j++) {
        size_t b = lzf.fact[j];
        size_t e = lzf.fact[j] + factLen(lzf, j);
        if (b <= ps2[i + missed].b - 1 && e >= ps2[i + missed].e - 1) {
          printf("%zu: ", j);
          fprintnf(stdout, lzf.str + lzf.fact[j], factLen(lzf, j));
          printf("\n");
          fact = lzf.fact[j];
        }
      }
      int64_t prev = fact;
      fact = 0;
      while (prev != -1) {
        fact = prev;
        prev = lzf.prevOcc[fact];
        printf("pos: %ld\n", prev);
      }

      missed++;
    }
  }

  mu_assert_eq(ps2.size(), ps.size(), "number of periodicities does not match");
  return NULL;
}

// compare: getting all, then filtering afterwards (as proposed in book)
// vs. getting only runs directly (for better performance)
char const *test_randomOnlyRuns() {
  size_t n = 1000000;
  string str = randSeq(n);
  char const *s = str.c_str();
  n++; //$ border
  fprintnf(stdout, s, 80);
  printf("\n");
  Esa esa(s, n); // calculate esa
  Fact lzf;
  computeLZFact(lzf, esa, false);

  size_t plen, plen2;
  // get all, then remove non-runs
  auto pl = getPeriodicityLists(false, lzf, esa, plen);
  for (size_t i = 0; i < lzf.strLen; i++) {
    Periodicity lastp(0, 0, 0);
    bool first = true;
    auto it = pl[i].begin();
    while (it != pl[i].end()) {
      Periodicity currp = *it;
      if (!first && lastp.b == currp.b && lastp.e == currp.e) {
        pl[i].erase(it++);
        plen--;
      } else {
        lastp = currp;
        it++;
      }
      first = false;
    }
  }

  // get just runs
  auto ps = getPeriodicities(true, lzf, esa, plen2);

  mu_assert_eq(plen2, plen, "number of runs does not match");
  return NULL;
}

char const *all_tests() {
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
