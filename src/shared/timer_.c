#include "timer_.h"

#include <time.h>
#include <sys/time.h>

long long get_thread_ns()
{
  struct timespec tp;
  if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tp))
    return -1LL;
  return (tp.tv_sec * 1000000000LL + tp.tv_nsec);
}

long long get_sys_us()
{
  struct timeval tv;
  if (gettimeofday(&tv, NULL))
    return -1LL;
  return (tv.tv_sec * 1000000LL + tv.tv_usec);
}
