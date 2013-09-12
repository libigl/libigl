#include "load_shader.h"
#ifndef IGL_NO_OPENGL

// Copyright Denis Kovacs 4/10/08
#include "print_shader_info_log.h"
#include <cstdio>
IGL_INLINE GLuint igl::load_shader(const char *src,const GLenum type)
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
