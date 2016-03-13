#include <sys/resource.h>
#include "args.h"

#define CPU_TIME                                                                         \
  (getrusage(RUSAGE_SELF, &ruse),                                                        \
   ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +                                         \
       1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))
void tick();
void tock(char *str);
