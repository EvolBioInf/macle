#pragma once
#include <cinttypes>
#include <cstdlib>
#include <vector>

#ifdef USE_SDSL
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wtype-limits"
#include <sdsl/int_vector.hpp>
#include <sdsl/rmq_support.hpp>
#pragma GCC diagnostic pop
#ifndef VECBIT
#define VECBIT 0
#endif
typedef sdsl::int_vector<VECBIT> uint_vec;
#else
typedef std::vector<uint32_t> uint_vec;
#endif
