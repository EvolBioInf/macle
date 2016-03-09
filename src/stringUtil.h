/***** stringUtil.h *********************************************************
 * Description:  Header file for string utilities.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Sun Jun  6 00:23:34 2004.
 * File created on Wed Dec  6 17:19:12 2006.
 *****************************************************************************/
#pragma once

char *randSeq(size_t n);

char *chomp(char *line);
char *cleanWord(char *word);
char *cleanWordEdges(char *word);

char *strdup2(char *s);

void strtolower(char *s, size_t l);
void strtoupper(char *s, size_t l);

void reverse(char *s);
void replace(char *string, char original, char replacement);
void split(char *line, char *splitC, char **splitArray, int *arrayLen);

void fprintnf(FILE *fp, char *str, int n);
