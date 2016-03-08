#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"
#include "interval.h"

#include <fcntl.h>
#include <unistd.h>

char *test_intervals() {
  int fd = open("Data/hotspotExample2.fasta", 0);
  Sequence *seq = readFasta(fd);
  close(fd);

  char *s = seqStr(seq,0);
  size_t n = seqLen(seq,0);
  Esa *esa = getEsa(s, n+1); //calculate esa, including $

  Interval iv;
  //do some perfect matches
  iv = getInterval(esa, "$", 1);
  mu_assert(iv.lb==0 && iv.rb==0 && iv.lcp==1, "wrong interval for $");

  iv = getInterval(esa, "AA", 2);
  mu_assert(iv.lb==2 && iv.rb==2 && iv.lcp==2, "wrong interval for AA");

  iv = getInterval(esa, "ACGC", 4);
  mu_assert(iv.lb==34 && iv.rb==35 && iv.lcp==4, "wrong interval for ACGC");

  iv = getInterval(esa, "CACA", 4);
  mu_assert(iv.lb==41 && iv.rb==71 && iv.lcp==4, "wrong interval for CACA");

  iv = getInterval(esa, "TGT", 3);
  mu_assert(iv.lb==98 && iv.rb==100 && iv.lcp==3, "wrong interval for TGT");

  //only a prefix matched
  iv = getInterval(esa, "GTGTA", 5);
  mu_assert(iv.lb==89 && iv.rb==90 && iv.lcp==4, "wrong interval for GTGTA");

  //nothing matched -> lcp 0, "complete" interval
  iv = getInterval(esa, "X", 1);
  mu_assert(iv.lb==0 && iv.rb==(int64_t)esa->n-1 && iv.lcp==0, "wrong interval for X");


  freeEsa(esa);
  freeSequence(seq);
  return NULL;
}

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_intervals);
  return NULL;
}
RUN_TESTS(all_tests)
