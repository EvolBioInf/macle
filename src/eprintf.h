/***** eprintf.h **************************************************
 * Description: Header file for eprintf, which provides error-
 *              handling capabilities.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * File created on Fri Dec 17 11:16:37 2004.
 *****************************************************************/
#pragma once

#include <stdio.h>

void setprogname2(char *);
void eprintf(char *, ...);

FILE *efopen(char *fname, char *mode);
int eopen(char *fname, int flag);
char *estrdup(char *);
void *emalloc(size_t);
void *ecalloc(size_t, size_t);
void *erealloc(void *, size_t);
