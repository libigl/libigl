#include "compile_and_link_program.h"
#include "compile_shader.h"
#include "report_gl_error.h"
#include <iostream>
#include <cassert>

IGL_INLINE GLuint igl::compile_and_link_program(
  const char * v_str, const char * f_str)
{
  GLuint vid = compile_shader(GL_VERTEX_SHADER,v_str);
  GLuint fid = compile_shader(GL_FRAGMENT_SHADER,f_str);

  GLuint prog_id = glCreateProgram();
  assert(prog_id != 0 && "Failed to create shader.");
  glAttachShader(prog_id,vid);
  igl::report_gl_error("glAttachShader (vid): ");
  glAttachShader(prog_id,fid);
  igl::report_gl_error("glAttachShader (fid): ");

  glLinkProgram(prog_id);
  igl::report_gl_error("glLinkProgram: ");

  GLint status;
  glGetProgramiv(prog_id, GL_LINK_STATUS, &status);
  if (status != GL_TRUE)
  {
    char buffer[512];
    glGetProgramInfoLog(prog_id, 512, NULL, buffer);
    std::cerr << "Linker error: " << std::endl << buffer << std::endl;
    prog_id = 0;
  }
  return prog_id;
}
