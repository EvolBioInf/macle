#include "prelude.h"
#include "eprintf.h"
#include "pfasta.h"
#include "fastafile.h"
#include <string.h>
#include <err.h>
#include <fcntl.h>
#include <unistd.h>

void free_fasta_file(FastaFile *ff) {
  for (size_t i = 0; i < ff->n; i++)
    pfasta_seq_free(&(ff->seq[i]));
  free(ff->seq);
  free(ff);
}

// input: filename of fasta file (or 0 for STDIN), reference to an empty number
// output: either successfully parsed file, or NULL
FastaFile *read_fasta_file(char *file) {
  // open file
  int fd = file ? eopen(file, O_RDONLY) : STDIN_FILENO;
  char *filename = file ? file : "stdin";
  if (fd < 0)
    err(1, "%s", file);

  // init parser
  int l;
  pfasta_file pf;
  if ((l = pfasta_parse(&pf, fd)) != 0) {
    warnx("%s: %s", filename, pfasta_strerror(&pf));
    return NULL;
  }

  // read sequences
  size_t n = 0;
  pfasta_seq *ps = NULL;
  pfasta_seq seq;
  while ((l = pfasta_read(&pf, &seq)) == 0) {
    ps = erealloc(ps, (++n) * sizeof(pfasta_seq));
    ps[n - 1] = seq;
  }

  if (l < 0) {
    pfasta_seq_free(&seq);
    warnx("%s: %s", filename, pfasta_strerror(&pf));
    for (size_t i = 0; i < n; i++)
      pfasta_seq_free(&ps[i]);
    free(ps);
    pfasta_free(&pf);
    return NULL;
  }

  pfasta_free(&pf);
  if (file)
    close(fd);

  FastaFile *ff = emalloc(sizeof(FastaFile));
  ff->seq = ps;
  ff->n = n;
  for (size_t i = 0; i < n; i++)
    ps[i].len = strlen(ps[i].seq) - 1;
  return ff;
}

/*
int main(int argc, char *argv[]) {
  FastaFile *ff = read_fasta_file(argv[1]);
  if (!ff) {
    printf("error while reading file\n");
    return 1;
  }

  printf("Read %ld sequences from %s:\n", ff->n, argv[1]);
  for (size_t i=0; i<ff->n; i++) {
    printf("%ld: %s -> %s\n", ff->seq[i].len, ff->seq[i].name, ff->seq[i].seq);
  }
}
*/
