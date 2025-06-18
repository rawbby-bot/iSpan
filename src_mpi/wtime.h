#ifndef __H_TIME__
#define __H_TIME__

#include <stdlib.h>
#include <sys/time.h>

inline double
wtime()
{
  double time[2];
  timeval time1{};
  gettimeofday(&time1, nullptr);

  time[0] = static_cast<double>(time1.tv_sec);
  time[1] = static_cast<double>(time1.tv_usec);

  return time[0] + time[1] * 1.0e-6;
}

#endif
