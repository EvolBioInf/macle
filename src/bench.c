#include "bench.h"

static struct rusage ruse;
static double last_tick_time = 0;

void tick() { last_tick_time = CPU_TIME; }

void tock(char *str) {
  if (args.b)
    fprintf(stderr, "[BENCH] %s: %.2fs\n", str, CPU_TIME - last_tick_time);
}
