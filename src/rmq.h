#pragma once
#include <cinttypes>
#include <cstdlib>
#include <vector>

#include <sdsl/int_vector.hpp>

class RMQ {
public:
  RMQ(sdsl::int_vector<VECBIT> const &A);
  int64_t get(size_t l, size_t r) const;

  sdsl::int_vector<VECBIT> const *arr;

  std::vector<int64_t> btab;
  std::vector<int64_t> ptab;
};

std::vector<int64_t> precomputePow2RMQ(std::vector<int64_t> const &A);
int64_t getRMQwithPow2(std::vector<int64_t> const &A, std::vector<int64_t> const &B,
                       size_t l, size_t r);
std::vector<int64_t> precomputeBlockRMQ(std::vector<int64_t> const &A);
