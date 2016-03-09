/***** stringUtil.c **********************************************
 * Description: Collection of string handling functions.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 10:02:16 2004.
 ****************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "stringUtil.h"
#include "eprintf.h"

/* generate random sequence of given length */
char *randSeq(size_t n) {
  static char *alphabet = "ACGT";
  char *s = emalloc((n + 2) * sizeof(char));
  s[n] = '$';
  s[n + 1] = '\0';
  for (size_t i = 0; i < n; i++)
    s[i] = alphabet[rand() % 4];
  return s;
}

/* chomp: remove carriage return from string */
char *chomp(char *line) {
  int i, l;
  l = strlen(line);
  for (i = 0; i < l; i++) {
    if (line[i] == '\n' || line[i] == '\r') {
      line[i] = '\0';
      break;
    }
  }
  return line;
}

/* fprintnf: print max of n characters of str onto fp; add ... if
 *   str was truncated
 */
void fprintnf(FILE *fp, char *str, int n) {
  int i, l, m;
  l = strlen(str);
  m = n < l ? n : l;
  for (i = 0; i < m; i++)
    fprintf(fp, "%c", str[i]);
  if (m < l)
    fprintf(fp, "...");
}

void strtolower(char *s, size_t l) {
  for (size_t i = 0; i < l; i++)
    s[i] = tolower(s[i]);
}

void strtoupper(char *s, size_t l) {
  for (size_t i = 0; i < l; i++)
    s[i] = toupper(s[i]);
}

/* strdup: make a duplicate of s */
char *strdup2(char *s) {
  char *p;
  p = (char *)malloc(strlen(s) + 1); /* +1 for '\0' */
  if (p != NULL)
    p = strcpy(p, s);
  return p;
}

void split(char *line, char *splitC, char *splitArray[], int *arrayLen) {
  char *line2, c;
  int i, j, start, end;

  start = end = 0;
  *arrayLen = 0;
  line2 = (char *)malloc(strlen(line) + 1); /* +1 for '\0' */
  if (line2 != NULL)
    line2 = strcpy(line2, line);
  i = 0;
  /* count number of fields */
  while ((c = line2[end++]) != '\0') {
    if (c == *splitC) {
      splitArray[*arrayLen] = (char *)malloc((end - start) * sizeof(char));
      i = 0;
      for (j = start; j < end - 1; j++) {
        splitArray[*arrayLen][i++] = line2[j];
      }
      splitArray[*arrayLen][i] = '\0';
      start = end;
      (*arrayLen)++;
    }
  }
  i = 0;
  splitArray[*arrayLen] = (char *)malloc((end - start + 1) * sizeof(char));
  for (j = start; j < end; j++) {
    splitArray[*arrayLen][i++] = line2[j];
  }
  splitArray[*arrayLen][i] = '\0';
  (*arrayLen)++;
}

void replace(char *string, char original, char replacement) {
  long i, l;
  l = strlen(string);
  for (i = 0; i < l; i++) {
    if (*string == original)
      *string = replacement;
    string++;
  }
}

/* remove leading and trailing non-alphanumerical characters */
char *cleanWordEdges(char *word) {
  int l1, l2;
  char *word2 = NULL;

  /* remove leading non-alphanumericals */
  while (!isalnum(*word) && *word != '\0') {
    word++;
  }

  /* remove trailing non-alphanumericals */
  l1 = strlen(word);
  l2 = l1;
  while (l2 > 0 && !isalnum(word[l2 - 1])) {
    l2--;
  }
  if (l2 < l1) {
    word2 = (char *)emalloc((l2 + 1) * sizeof(char));
    word2 = strncpy(word2, word, l2);
    return word2;
  } else {
    return word;
  }
}

/* remove all non-alphanumerical characters */
char *cleanWord(char *word) {
  int i, j, l1, l2;
  char *word2 = NULL;

  /* get word length */
  l1 = strlen(word);

  /* count alphanumericals */
  l2 = 0;
  for (i = 0; i < l1; i++)
    if (isalnum(word[i]))
      l2++;
  if (l2 < l1) { /* create cleaned string */
    word2 = (char *)emalloc((l2 + 1) * sizeof(char));
    j = 0;
    for (i = 0; i < l2; i++)
      if (isalnum(word[i]))
        word2[j++] = word[i];
    word2[j] = '\0';
    return word2;
  } else {
    return word;
  }
}

/* reverse: reverse stirng s in place */
void reverse(char *s) {
  int c, i, j;
  for (i = 0, j = strlen(s) - 1; i < j; i++, j--) {
    c = s[i];
    s[i] = s[j];
    s[j] = c;
  }
}
