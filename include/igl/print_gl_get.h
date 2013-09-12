#ifndef IGL_PRINT_GL_GET_H
#define IGL_PRINT_GL_GET_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

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
#endif
