#pragma once
#include <cinttypes>
#include <cstdlib>

class RMQ {
public:
  RMQ(int64_t *A, size_t n);
  ~RMQ();
  int64_t get(size_t l, size_t r) const;

  int64_t *arr;
  const size_t n;

  int64_t *tab;
};

void precomputePow2RMQ(int64_t *A, size_t n, int64_t *B);
int64_t getRMQwithPow2(int64_t *A, size_t n, int64_t *B, size_t l, size_t r);
void precomputeBlockRMQ(int64_t *A, size_t n, int64_t *B);

int64_t RMQuery(int64_t const *A, size_t n, int64_t const *B, size_t l, size_t r);
