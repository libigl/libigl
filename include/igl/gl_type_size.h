#ifndef IGL_GL_TYPE_SIZE_H
#define IGL_GL_TYPE_SIZE_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

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
#endif
