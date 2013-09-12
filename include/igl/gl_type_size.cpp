#include "gl_type_size.h"
#ifndef IGL_NO_OPENGL
#include <cassert>

IGL_INLINE int igl::gl_type_size(const GLenum type)
{
  switch(type)
  {
    case GL_DOUBLE:
      return 8;
      break;
    case GL_FLOAT:
      return 4;
      break;
    case GL_INT:
      return 4;
      break;
    default:
      // should handle all other GL_[types]
      assert(false);
      break;
  }
  return -1;
}
#endif
