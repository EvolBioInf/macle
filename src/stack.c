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

Stack *newStack(size_t max) {
  Stack *s = emalloc(sizeof(Stack));
  s->currMax = max;
  s->array = emalloc(max * sizeof(stackel));
  s->n = 0;
  return s;
}

void freeStack(Stack *s) {
  free(s->array);
  free(s);
}

bool stackEmpty(Stack *s) { return s->n == 0; }

void stackPush(Stack *s, stackel stackItem) {
  if (s->n == s->currMax) {
    s->currMax *= 2;
    s->array = (stackel *)erealloc(s->array, s->currMax * sizeof(stackel));
  }
  s->array[s->n++] = stackItem;
}

stackel stackTop(Stack *s) {
  if (s->n > 0)
    return s->array[s->n - 1];
  else {
    printf("ERROR[stack]: trying to pop an empty stack.\n");
    return (void *)-1;
  }
}

stackel stackPop(Stack *s) {
  if (s->n > 0)
    return s->array[--(s->n)];
  else {
    printf("ERROR[stack]: trying to pop an empty stack.\n");
    return (void *)-1;
  }
}
