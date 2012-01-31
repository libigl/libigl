#ifndef IGL_LOAD_SHADER_H 
#define IGL_LOAD_SHADER_H
#include "igl_inline.h" 

#ifdef __APPLE__
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
  // Creates and compiles a shader from a given string
  // Inputs:
  //   src  string containing GLSL shader code
  //   type  GLSL type of shader, one of:
  //     GL_VERTEX_SHADER
  //     GL_FRAGMENT_SHADER
  //     GL_GEOMETRY_SHADER
  // Returns  index id of the newly created shader, 0 on error
  IGL_INLINE GLuint load_shader(const char *src,const GLenum type);
}

#ifdef IGL_HEADER_ONLY
#  include "load_shader.cpp"
#endif

#endif
