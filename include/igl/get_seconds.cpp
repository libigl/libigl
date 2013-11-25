// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "get_seconds.h"
// NULL for Linux
#include <cstddef>

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
