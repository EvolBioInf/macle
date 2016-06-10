#include "fastafile.h"
#include "util.h"

#include <err.h>
#include "pfasta.h"

//for open_or_fail flags
#include <unistd.h>
#include <fcntl.h>

#include <cstring>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

FastaSeq::FastaSeq() {}
FastaSeq::FastaSeq(string const &n, string const &c, string const &s)
    : name(n), comment(c), seq(s) {}

FastaFile::FastaFile() : failed(false) {}

// input: filename of fasta file (or 0 for STDIN), reference to an empty number
// output: either successfully parsed file, or NULL
FastaFile::FastaFile(char const *file) : failed(false) {
  seqs = vector<FastaSeq>();

  // open file
  int fd = file ? open_or_fail(file, O_RDONLY) : STDIN_FILENO;
  filename = file ? base_name(string(file)) : "<stdin>";
  if (fd < 0)
    err(1, "%s", file);

  // init parser
  int l;
  pfasta_file pf;
  if ((l = pfasta_parse(&pf, fd)) != 0) {
    warnx("%s: %s", filename.c_str(), pfasta_strerror(&pf));
    pfasta_free(&pf);
    failed = true;
    return;
  }

  // read sequences
  pfasta_seq seq;
  while ((l = pfasta_read(&pf, &seq)) == 0) {
    seqs.push_back(FastaSeq());
    auto &fseq = seqs.back();
    fseq.name = seq.name ? string(seq.name) : "";
    fseq.comment = seq.comment ? string(seq.comment) : "";
    fseq.seq = seq.seq ? string(seq.seq) : "";
    for (char &c : fseq.seq)
      c = toupper(c); // acgt->ACGT
    pfasta_seq_free(&seq);
  }

  if (l < 0) {
    warnx("%s: %s", filename.c_str(), pfasta_strerror(&pf));
    pfasta_seq_free(&seq);
    failed = true;
  }

  pfasta_free(&pf);
  if (file)
    close(fd);
}
