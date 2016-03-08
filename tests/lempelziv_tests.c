#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"
#include "lempelziv.h"

#include <fcntl.h>
#include <unistd.h>

char *factors[] = {
 "G","C","A","C","GCACGCAC",
 "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
 "T","AT","GC","TA","AC","T","CTC","A","G","TCT","GT","GTGTG","CA","$"
};

char *test_LempelZiv() {
  int fd = open("Data/hotspotExample2.fasta", 0);
  Sequence *seq = readFasta(fd);
  close(fd);

  char *s = getSeq(seq,0);
  size_t n = seqLen(seq,0);
  Esa *esa = getEsa(s, n+1); //calculate esa, including $

  Fact *lzf = computeLZFact(esa);
  mu_assert(lzf->n == 20, "wrong number of LZ factors");

  for (size_t i=0; i<lzf->n; i++) {
    mu_assert(!strncmp(lzf->str+lzf->fact[i], factors[i], factLen(lzf, i)), "wrong factor");
  }

  freeEsa(esa);
  freeSequence(seq);
  return NULL;
}

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_LempelZiv);
  return NULL;
}
RUN_TESTS(all_tests)
