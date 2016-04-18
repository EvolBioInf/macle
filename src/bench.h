#include <sys/resource.h>
#include "args.h"

#define CPU_TIME                                                                         \
  (getrusage(RUSAGE_SELF, &ruse),                                                        \
   ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +                                         \
       1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

static struct rusage ruse;
static double last_tick_time = 0;

#define tick()                                                                           \
  do {                                                                                   \
    last_tick_time = CPU_TIME;                                                           \
  } while (0);

#define tock(str)                                                                        \
  do {                                                                                   \
    if (args.b)                                                                          \
      fprintf(stderr, "[BENCH] %s: %.2fs\n", str, CPU_TIME - last_tick_time);            \
  } while (0);
