#include "minunit.h"

#include "stack.h"
#include "list.h"

char *test_list() {
  List *l = NULL;
  for (int64_t i = 500; i < 1000; i++)
    listAppend(&l, (void *)i);
  for (int64_t i = 499; i >= 0; i--)
    listPrepend(&l, (void *)i);

  List *curr = l;
  size_t i = 0;
  do {
    mu_assert(curr->value == (void *)i, "wrong list value");
    i++;
    curr = curr->next;
  } while (curr);

  freeList(&l);
  mu_assert(l == NULL, "list pointer not NULL after free");

  return NULL;
}

char *test_stack() {
  Stack *s = newStack(10);
  mu_assert(stackEmpty(s), "new stack not empty!");
  mu_assert(stackTop(s) == (void *)-1, "top of empty stack wrong");
  mu_assert(stackPop(s) == (void *)-1, "pop of empty stack wrong");
  for (int64_t i = 0; i < 1000; i++)
    stackPush(s, (stackel)i);
  mu_assert(s->n == 1000, "stack size wrong");
  mu_assert((int64_t)stackTop(s) == 999, "top element wrong");
  mu_assert((int64_t)stackTop(s) == 999, "top element changed after stackTop");

  for (int64_t i = 999; i >= 0; i--)
    mu_assert((int64_t)stackPop(s) == i, "pop of stack wrong");
  mu_assert(stackEmpty(s), "stack still not empty");

  freeStack(s);
  return NULL;
}

char *all_tests() {
  mu_suite_start();

  mu_run_test(test_list);
  mu_run_test(test_stack);

  return NULL;
}

RUN_TESTS(all_tests)
