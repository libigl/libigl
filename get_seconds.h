#include <sys/time.h>

namespace igl
{
  // Return the current time since Epoch in seconds
  double get_seconds();
}

//Implementation

double igl::get_seconds()
{
  timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec + time.tv_usec / 1e6;
}
