#pragma once
#include <cinttypes>
#include <cstdlib>
#include <vector>

#ifdef USE_SDSL
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#ifndef VECBIT
#define VECBIT 0
#endif
typedef sdsl::int_vector<VECBIT> uint_vec;
#else
typedef std::vector<uint32_t> uint_vec;
#endif

class RMQ_custom {
public:
  RMQ_custom(uint_vec const *A);
  size_t operator()(size_t l, size_t r) const;

  uint_vec const *arr;

  uint_vec btab;
  uint_vec ptab;
};

uint_vec precomputePow2RMQ(uint_vec const &A);
int64_t getRMQwithPow2(std::vector<int64_t> const &A, std::vector<int64_t> const &B,
                       size_t l, size_t r);
uint_vec precomputeBlockRMQ(uint_vec const &A);

#ifdef USE_SDSL
typedef sdsl::rmq_succinct_sct<> RMQ;
#else
typedef RMQ_custom RMQ;
#endif
