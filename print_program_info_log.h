#ifndef IGL_PRINT_PROGRAM_INFO_LOG_H
#define IGL_PRINT_PROGRAM_INFO_LOG_H
#include "igl_inline.h"

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#endif
namespace igl
{
  // Inputs:
  //   obj  OpenGL index of program to print info log about
  IGL_INLINE void print_program_info_log(const GLuint obj);
}

#ifdef IGL_HEADER_ONLY
#  include "print_program_info_log.cpp"
#endif

#endif
