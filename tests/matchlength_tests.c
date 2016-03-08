#include "minunit.h"

#include "sequenceData.h"
#include "esa.h"
#include "matchlength.h"

#include <fcntl.h>
#include <unistd.h>

static char *factors1[] = {"CCC","C","G","CTC","TC","C","A","$"};
static char *factors2[] = { "GCACGCAC","GCAC",
  "ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA",
  "TA","TGC","TA","AC","TCT","CA","GT","CT","GTGTG","TGC","A","$"
};

char *checkML(char *file, char *facts[], size_t num) {
  int fd = open(file, 0);
  Sequence *seq = readFasta(fd);
  close(fd);

  char *s = seqStr(seq,0);
  size_t n = seqLen(seq,0);
  Esa *esa = getEsa(s, n+1); //calculate esa, including $

  Fact *mlf = computeMLFact(esa);
  mu_assert(mlf->n == num, "wrong number of ML factors");

  for (size_t i=0; i<mlf->n; i++) {
    mu_assert(!strncmp(mlf->str+mlf->fact[i], facts[i], factLen(mlf, i)), "wrong factor");
  }

  freeFact(mlf);
  freeEsa(esa);
  freeSequence(seq);
  return NULL;
}

char *test_MatchLength1() { return checkML("Data/hotspotExample1.fasta", factors1, 8); }
char *test_MatchLength2() { return checkML("Data/hotspotExample2.fasta", factors2, 15); }

char *all_tests() {
  mu_suite_start();
  mu_run_test(test_MatchLength1);
  mu_run_test(test_MatchLength2);
  return NULL;
}
RUN_TESTS(all_tests)
