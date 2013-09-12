#ifndef IGL_UNIFORM_TYPE_TO_STRING_H
#define IGL_UNIFORM_TYPE_TO_STRING_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include <string>

#include "OpenGL_convenience.h"

namespace igl
{
  // Convert a GL uniform variable type (say, returned from
  // glGetActiveUniform) and output a string naming that type
  // Inputs:
  //   type  enum for given type
  // Returns string name of that type
  IGL_INLINE std::string uniform_type_to_string(const GLenum type);
}

#ifdef IGL_HEADER_ONLY
#  include "uniform_type_to_string.cpp"
#endif

#endif
#endif
