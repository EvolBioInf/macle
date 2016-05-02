#include <ctime>
#include "minunit.h"
#include "rmq.h"

#include <vector>

char const *test_rmqSmall() {
  size_t n = 25;
  std::vector<int64_t> array(n);
  for (size_t i = 0; i < n; i++) {
    array[i] = rand() % n;
    /* printf("%ld ", array[i]); */
  }
  /* printf("\n\n"); */

  RMQ rmq(array);

  for (size_t i = 0; i < n; i++)
    for (size_t j = i; j < n; j++) {
      int64_t exp = INT64_MAX;
      for (size_t k = i; k <= j; k++)
        if (array[k] < exp)
          exp = array[k];

      int64_t obs = rmq.get(i, j);
      if (exp != obs)
        printf("tried %zu and %zu\n", i, j);
      mu_assert_eq(exp, obs, "wrong range minimum");
    }

  return NULL;
}

char const *test_rmqRand() {
  size_t n = 10000;
  std::vector<int64_t> array(n);
  for (size_t i = 0; i < n; i++) {
    array[i] = rand() % n;
  }

  RMQ rmq(array);

  for (size_t i = 0; i < 1000; i++) {
    size_t l = rand() % n;
    size_t r = l + rand() % (n - l);
    int64_t exp = INT64_MAX;
    for (size_t k = l; k <= r; k++)
      if (array[k] < exp)
        exp = array[k];

    int64_t obs = rmq.get(l, r);
    if (exp != obs)
      printf("tried %zu and %zu\n", l, r);
    mu_assert_eq(exp, obs, "wrong range minimum");
  }

  return NULL;
}

char const *all_tests() {
  srand(time(NULL));
  mu_suite_start();
  mu_run_test(test_rmqSmall);
  for (size_t i = 0; i < 3; i++)
    mu_run_test(test_rmqRand);
  return NULL;
}
RUN_TESTS(all_tests)
