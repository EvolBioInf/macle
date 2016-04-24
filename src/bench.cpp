#include <chrono>
#include <cstdio>
using namespace std::chrono;

#include "args.h"

static high_resolution_clock::time_point tp;

void tick() { tp = high_resolution_clock::now(); }

void tock(char const *str) {
  if (args.b) {
    auto const tp2 = high_resolution_clock::now();
    double span = duration_cast<duration<double>>(tp2 - tp).count();
    fprintf(stderr, "[BENCH] %s: %.2fs\n", str, span);
  }
}
