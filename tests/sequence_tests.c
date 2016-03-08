#include "minunit.h"

#include "sequenceData.h"
#include <string.h>

#include <fcntl.h>
#include <unistd.h>

char *test_readFasta() {
  int fd = open("Data/seqTest01.fasta", 0);
  Sequence *seq = readFasta(fd);
  close(fd);

  mu_assert(seq->numSeq == 5, "wrong number of sequences identified");
  mu_assert(seqLen(seq, 0) == 1, "wrong sequence length");
  mu_assert(seqLen(seq, 1) == 4, "wrong sequence length");
  mu_assert(seqLen(seq, 2) == 6, "wrong sequence length");
  mu_assert(seqLen(seq, 3) == 0, "wrong sequence length");
  mu_assert(seqLen(seq, 4) == 3, "wrong sequence length");

  mu_assert(!strcmp(seq->headers[4], ">SeqTest01e"), "wrong sequence name");
  mu_assert(getSeq(seq,2)[seqLen(seq,2)] == '$', "border not present");

  freeSequence(seq);
  return NULL;
}

char *all_tests() {
  mu_suite_start();

  mu_run_test(test_readFasta);

  return NULL;
}

RUN_TESTS(all_tests)
