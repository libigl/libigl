#ifndef IGL_GET_SECONDS_H
#define IGL_GET_SECONDS_H

namespace igl
{
  // Return the current time in seconds since program start
  inline double get_seconds();

}

//Implementation
#if _WIN32
#  include <ctime>
inline double igl::get_seconds()
{
  // This does not work on mac os x with glut in the main loop
  return double(clock())/CLOCKS_PER_SEC;
}
#else
#  include <sys/time.h>
inline double igl::get_seconds()
{
  timeval time;
  gettimeofday(&time, NULL);
  return time.tv_sec + time.tv_usec / 1e6;
  // This does not work on mac os x with glut in the main loop
  //return double(clock())/CLOCKS_PER_SEC;
}
#endif
#endif
