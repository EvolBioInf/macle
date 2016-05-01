#include "fastafile.h"

#include <err.h>
#include "pfasta.h"

#include <errno.h>
#include <unistd.h>
#include <fcntl.h>

#include <cstring>
#include <iostream>
#include <vector>
using namespace std;

FastaFile::~FastaFile() {
  if (failed)
    return;
  for (auto it = seqs.begin(); it != seqs.end(); it++)
    pfasta_seq_free(&(*it));
}

/* eopen: open file on system level and report on error */
int open_or_fail(char const *fname, int flag) {
  int fd = open(fname, flag, 0);
  if (fd < 0) {
    cerr << "open(" << fname << "," << flag << ") failed: " << strerror(errno) << endl;
    exit(EXIT_FAILURE);
  }
  return fd;
}

// input: filename of fasta file (or 0 for STDIN), reference to an empty number
// output: either successfully parsed file, or NULL
FastaFile::FastaFile(char const *file) : failed(false) {
  seqs = vector<pfasta_seq>();

  // open file
  int fd = file ? open_or_fail(file, O_RDONLY) : STDIN_FILENO;
  char const *filename = file ? file : "stdin";
  if (fd < 0)
    err(1, "%s", file);

  // init parser
  int l;
  pfasta_file pf;
  if ((l = pfasta_parse(&pf, fd)) != 0) {
    warnx("%s: %s", filename, pfasta_strerror(&pf));
    this->failed = true;
    pfasta_free(&pf);
    return;
  }

  // read sequences
  pfasta_seq seq;
  while ((l = pfasta_read(&pf, &seq)) == 0)
    seqs.push_back(seq);

  if (l < 0) {
    pfasta_seq_free(&seq);
    warnx("%s: %s", filename, pfasta_strerror(&pf));
    for (auto it = seqs.begin(); it != seqs.end(); it++)
      pfasta_seq_free(&(*it));
    seqs.clear();
    pfasta_free(&pf);
    this->failed = true;
    return;
  }

  pfasta_free(&pf);
  if (file)
    close(fd);

  for (auto it = seqs.begin(); it != seqs.end(); it++)
    it->len = strlen(it->seq) - 1;
}
