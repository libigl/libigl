#ifndef IGL_GET_SECONDS_H
#define IGL_GET_SECONDS_H

namespace igl
{
  // Return the current time since Epoch in seconds
  inline double get_seconds();

}

//Implementation
#include <sys/time.h>

inline double igl::get_seconds()
{
  timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec + time.tv_usec / 1e6;
}
#endif
