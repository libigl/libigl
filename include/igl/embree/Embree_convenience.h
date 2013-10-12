#ifndef IGL_EMBREE_CONVENIENCE_H
#define IGL_EMBREE_CONVENIENCE_H

#undef interface
#undef near
#undef far
// Why are these in quotes? isn't that a bad idea?
#ifdef __GNUC__
// This is how it should be done
#  if __GNUC__ >= 4
#    if __GNUC_MINOR__ >= 6
#      pragma GCC diagnostic push
#      pragma GCC diagnostic ignored "-Weffc++"
#    endif
#  endif
// This is a hack
#  pragma GCC system_header
#endif
#include <common/intersector.h>
#include <common/accel.h>
#ifdef __GNUC__
#  if __GNUC__ >= 4
#    if __GNUC_MINOR__ >= 6
#      pragma GCC diagnostic pop
#    endif
#  endif
#endif

#endif
