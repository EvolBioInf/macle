/***** eprintf.c **************************************************
 * Description: Collection of functions for error handling.
 * Reference: Kernighan, B. W. and Pike, R. (1999). The Practice
 *            of programming. Addision Wesley; chapter 4.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Dec 17 11:16:34 2004.
 *****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include "eprintf.h"

char const *progname = NULL; /* program name for messages */

/* eprintf: print error message and exit */
void eprintf(char const *fmt, ...) {
  va_list args;
  fflush(stdout);
  if (progname != NULL)
    fprintf(stderr, "%s: ", progname);

  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);

  if (fmt[0] != '\0' && fmt[strlen(fmt) - 1] == ':')
    fprintf(stderr, " %s", strerror(errno));
  fprintf(stderr, "\n");
  exit(EXIT_FAILURE);
}

/* efopen: open file and report if error */
FILE *efopen(char const *fname, char const *mode) {
  FILE *fp = fopen(fname, mode);
  if (fp == NULL)
    eprintf("efopen(%s, %s) failed:", fname, mode);
  return fp;
}

/* eopen: open file on system level and report if error */
int eopen(char const *fname, int flag) {
  int fd = open(fname, flag, 0);
  if (fd < 0)
    eprintf("eopen(%s, %d) failed:", fname, flag);
  return fd;
}

/* estrdup: duplicate a string, report if error */
char *estrdup(char *s) {
  char *t = (char *)malloc(strlen(s) + 1);
  if (t == NULL)
    eprintf("estrdup(\"%.20s\") failed:", s);
  strcpy(t, s);
  return t;
}

/* emalloc: malloc and report if error */
void *emalloc(size_t n) {
  void *p = malloc(n);
  if (p == NULL)
    eprintf("malloc of %u bytes failed:", n);
  return p;
}

/* ecalloc: calloc and report if error */
void *ecalloc(size_t n, size_t sz) {
  void *p = calloc(n, sz);
  if (p == NULL)
    eprintf("calloc of %u bytes failed:", n*sz);
  return p;
}

/* erealloc: realloc and report if error */
void *erealloc(void *p, size_t n) {
  p = realloc(p, n);
  if (p == NULL)
    eprintf("realloc of %u bytes failed:", n);
  return p;
}
