#ifndef IGL_GL_TYPE_SIZE_H
#define IGL_GL_TYPE_SIZE_H
#include "igl_inline.h"

#if __APPLE__
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
  // Return the number of bytes for a given OpenGL type
  // Inputs:
  //   type  enum value of opengl type
  // Returns size in bytes of type
  IGL_INLINE int gl_type_size(const GLenum type);
}

#ifdef IGL_HEADER_ONLY
#  include "gl_type_size.cpp"
#endif

#endif
