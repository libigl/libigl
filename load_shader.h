#ifndef IGL_LOAD_SHADER_H 
#define IGL_LOAD_SHADER_H 

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#else
#   include <GL/gl.h>
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
  GLuint load_shader(const char *src,const GLenum type);
}

// Implmentation
// Copyright Denis Kovacs 4/10/08
#include "print_shader_info_log.h"
#include <cstdio>
GLuint igl::load_shader(const char *src,const GLenum type)
{
  GLuint s = glCreateShader(type);
  if(s == 0)
  {
    fprintf(stderr,"Error: load_shader() failed to create shader.\n");
    return 0;
  }
  // Pass shader source string
  glShaderSource(s, 1, &src, NULL);
  glCompileShader(s);
  // Print info log (if any)
  igl::print_shader_info_log(s);
  return s;
}

#endif

