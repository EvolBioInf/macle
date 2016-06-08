#include <ctime>
#include "minunit.h"
#include "rmq.h"

#include <sdsl/int_vector.hpp>
using namespace sdsl;

void test_rmqSmall() {
  size_t n = 25;
  int_vector<VECBIT> array(n);
  for (size_t i = 0; i < n; i++) {
    array[i] = rand() % n;
    /* printf("%ld ", array[i]); */
  }
  /* printf("\n\n"); */

  RMQ rmq(array);
  // rmq_succinct_sct<> rmq(&array);

  for (size_t i = 0; i < n; i++)
    for (size_t j = i; j < n; j++) {
      int64_t exp = INT64_MAX;
      for (size_t k = i; k <= j; k++)
        if ((int64_t)array[k] < exp) {
          exp = array[k];
        }

      // int64_t obs = array[rmq(i, j)];
      int64_t obs = rmq(i, j);
      mu_assert_eq(exp, obs, "wrong range minimum (tried " << i << " and " << j << ")");
    }
}

void test_rmqRand() {
  size_t n = 10000;
  sdsl::int_vector<VECBIT> array(n);
  for (size_t i = 0; i < n; i++) {
    array[i] = rand() % n;
  }

  RMQ rmq(array);
  // rmq_succinct_sct<> rmq(&array);

  for (size_t i = 0; i < 1000; i++) {
    size_t l = rand() % n;
    size_t r = l + rand() % (n - l);
    int64_t exp = INT64_MAX;
    for (size_t k = l; k <= r; k++)
      if ((int64_t)array[k] < exp)
        exp = array[k];

    // int64_t obs = array[rmq(l, r)];
    int64_t obs = rmq(l, r);
    mu_assert_eq(exp, obs, "wrong range minimum (tried " << l << " and " << r << ")");
  }
}

void all_tests() {
  srand(time(NULL));
  mu_run_test(test_rmqSmall);
  for (size_t i = 0; i < 3; i++)
    mu_run_test(test_rmqRand);
}
RUN_TESTS(all_tests)
