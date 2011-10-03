#ifndef IGL_GL_TYPE_SIZE
#define IGL_GL_TYPE_SIZE

#if __APPLE__
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

namespace igl
{
  // Return the number of bytes for a given OpenGL type
  // Inputs:
  //   type  enum value of opengl type
  // Returns size in bytes of type
  inline int gl_type_size(GLenum type);
}

//Implementation
inline int igl::gl_type_size(GLenum type)
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

