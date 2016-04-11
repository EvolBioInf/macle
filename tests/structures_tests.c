#include "minunit.h"

#include "kvec.h"
#include "list.h"

char *test_list() {
  List *l = NULL;
  for (int64_t i = 500; i < 1000; i++)
    listAppend(&l, (void *)i);
  for (int64_t i = 499; i >= 0; i--)
    listPrepend(&l, (void *)i);

  size_t i = 0;
  for (eachListItem(curr, l)) {
    mu_assert(curr->value == (void *)i, "wrong list value");
    i++;
  }

  mu_assert_eq(1000, listLength(l), "wrong list length");

  freeList(&l);
  mu_assert(l == NULL, "list pointer not NULL after free");

  return NULL;
}

char *test_stack() {
  kvec_t(int64_t) s;
  kv_init(s);
  mu_assert(kv_empty(s), "new stack not empty!");
  for (int64_t i = 0; i < 1000; i++)
    kv_push(int64_t, s, i);
  mu_assert(kv_size(s) == 1000, "stack size wrong");
  mu_assert(kv_top(s) == 999, "top element wrong");
  mu_assert(kv_top(s) == 999, "top element changed after stackTop");

  for (int64_t i = 999; i >= 0; i--)
    mu_assert(kv_pop(s) == i, "pop of stack wrong");
  mu_assert(kv_empty(s), "stack still not empty");

  kv_destroy(s);
  return NULL;
}

char *all_tests() {
  mu_suite_start();

  mu_run_test(test_list);
  mu_run_test(test_stack);

  return NULL;
}

RUN_TESTS(all_tests)
