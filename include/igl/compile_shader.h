#ifndef IGL_COMPILE_SHADER_H
#define IGL_COMPILE_SHADER_H
#include "OpenGL_convenience.h"
#include "igl_inline.h"
namespace igl
{
  // Compile a shader given type and string of shader code
  // 
  // Inputs:
  //   type  either GL_VERTEX_SHADER or GL_FRAGMENT_SHADER
  //   str  contents of shader code
  // Returns result of glCreateShader (id of shader)
  //
  // Example:
  //     GLuint vid = compile_shader(GL_VERTEX_SHADER,vertex_shader.c_str());
  //     GLuint fid = compile_shader(GL_FRAGMENT_SHADER,fragment_shader.c_str());
  //     GLuint prog_id = glCreateProgram();
  //     glAttachShader(prog_id,vid);
  //     glAttachShader(prog_id,fid);
  //     glLinkProgram(prog_id);
  //
  // Known bugs: seems to be duplicate of `load_shader`
  IGL_INLINE GLuint compile_shader(const GLint type, const char * str);
}
#ifndef IGL_STATIC_LIBRARY
#  include "compile_shader.cpp"
#endif
#endif
