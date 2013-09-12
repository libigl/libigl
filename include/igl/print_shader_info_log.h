#ifndef IGL_PRINT_SHADER_INFO_LOG_H
#define IGL_PRINT_SHADER_INFO_LOG_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

namespace igl
{
  // Inputs:
  //   obj  OpenGL index of shader to print info log about
  IGL_INLINE void print_shader_info_log(const GLuint obj);
}

#ifdef IGL_HEADER_ONLY
#  include "print_shader_info_log.cpp"
#endif

#endif
#endif
