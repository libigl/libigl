#ifndef IGL_UNIFORM_TYPE_TO_STRING_H
#define IGL_UNIFORM_TYPE_TO_STRING_H
#include "igl_inline.h"

#include <string>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#elif defined(_WIN32)
#  define NOMINMAX
#  include <Windows.h>
#  undef NOMINMAX
#  include <GL/glew.h>
#  include <GL/gl.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#endif

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
