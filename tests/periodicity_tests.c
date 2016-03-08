
#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"
#include "lempelziv.h"
#include "periodicity.h"

#include <fcntl.h>
#include <unistd.h>

char *test_periodicity() {
  int fd = open("Data/hotspotExample2.fasta", 0);
  Sequence *seq = readFasta(fd);
  close(fd);

  char *s = seqStr(seq,0);
  size_t n = seqLen(seq,0);
  Esa *esa = getEsa(s, n+1); //calculate esa, including $
  Fact *lzf = computeLZFact(esa);

  /* mu_assert(1==0, "TODO"); */

  freeFact(lzf);
  freeEsa(esa);
  freeSequence(seq);
  return NULL;
}

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_periodicity);
  return NULL;
}
RUN_TESTS(all_tests)
