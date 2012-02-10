#include "get_seconds.h"

#if _WIN32
#  include <ctime>
IGL_INLINE double igl::get_seconds()
{
  // This does not work on mac os x with glut in the main loop
  return double(clock())/CLOCKS_PER_SEC;
}
#else
#  include <sys/time.h>
IGL_INLINE double igl::get_seconds()
{
  timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec + time.tv_usec / 1e6;
  // This does not work on mac os x with glut in the main loop
  //return double(clock())/CLOCKS_PER_SEC;
}
#endif
