#ifndef IGL_PRINT_SHADER_INFO_LOG_H
#define IGL_PRINT_SHADER_INFO_LOG_H
#include "igl_inline.h"

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#elif defined(_WIN32)
#  define NOMINMAX
#  include <Windows.h>
#  undef NOMINMAX
#  include <GL/gl.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#endif

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
