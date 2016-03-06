/***** stack.c ************************************
 * Description: Array implementation of pushdown
 *   stack. Taken from p. 146 of Sedgewick, R.
 *   (1998). Algorithms in C. Third Edition.
 *   Addision-Wesely, Boston.
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Wed Nov 25 17:32:47 2009
 **************************************************/
#include "prelude.h"

#include "eprintf.h"
#include "stack.h"

static int *array;
static size_t n;
static size_t maxN;

void stackInit(size_t m) {
  maxN = m;
  array = (int *)emalloc(maxN * sizeof(int));
  n = 0;
}

bool stackEmpty() { return n == 0; }

void stackPush(int stackItem) {
  if (n == maxN) {
    maxN *= 2;
    array = (int *)erealloc(array, maxN * sizeof(int));
  }
  array[n++] = stackItem;
}

int stackTop() {
  if (n > 0)
    return array[n - 1];
  else {
    printf("ERROR[stack]: trying to pop an empty stack.\n");
    return -1;
  }
}

int stackPop() {
  if (n > 0)
    return array[--n];
  else {
    printf("ERROR[stack]: trying to pop an empty stack.\n");
    return -1;
  }
}

void freeStack() { free(array); }
