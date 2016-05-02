#pragma once
#include <cinttypes>
#include <cstdlib>
#include <vector>

class RMQ {
public:
  RMQ(std::vector<int64_t> const &A);
  int64_t get(size_t l, size_t r) const;

  std::vector<int64_t> const *arr;

  std::vector<int64_t> btab;
  std::vector<int64_t> ptab;
};

std::vector<int64_t> precomputePow2RMQ(std::vector<int64_t> const &A);
int64_t getRMQwithPow2(std::vector<int64_t> const &A, std::vector<int64_t> const &B, size_t l, size_t r);
std::vector<int64_t> precomputeBlockRMQ(std::vector<int64_t> const &A);
