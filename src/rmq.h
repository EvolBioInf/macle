#pragma once
#include <cinttypes>
#include <cstdlib>
#include <vector>

#include <sdsl/int_vector.hpp>

class RMQ {
public:
  RMQ(sdsl::int_vector<VECBIT> const &A);
  size_t operator()(size_t l, size_t r) const;

  sdsl::int_vector<VECBIT> const *arr;

  sdsl::int_vector<VECBIT> btab;
  sdsl::int_vector<VECBIT> ptab;
};

sdsl::int_vector<VECBIT> precomputePow2RMQ(std::vector<int64_t> const &A);
int64_t getRMQwithPow2(std::vector<int64_t> const &A, std::vector<int64_t> const &B,
                       size_t l, size_t r);
sdsl::int_vector<VECBIT> precomputeBlockRMQ(std::vector<int64_t> const &A);
