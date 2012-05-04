#ifndef IGL_PRINT_GL_GET_H
#define IGL_PRINT_GL_GET_H
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
  // Prints the value of pname found by issuing glGet*(pname,*)
  // Inputs:
  //   pname  enum key to gl parameter
  IGL_INLINE void print_gl_get(GLenum pname);
}

#ifdef IGL_HEADER_ONLY
#  include "print_gl_get.cpp"
#endif

#endif
