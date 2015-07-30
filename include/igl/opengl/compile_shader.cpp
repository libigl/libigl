// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "compile_shader.h"
#include "report_gl_error.h"
#include <iostream>

IGL_INLINE GLuint igl::opengl::compile_shader(const GLint type, const char * str)
{
  GLuint id = glCreateShader(type);
  report_gl_error("glCreateShader: ");
  glShaderSource(id,1,&str,NULL);
  report_gl_error("glShaderSource: ");
  glCompileShader(id);
  report_gl_error("glCompileShader: ");

  GLint status;
  glGetShaderiv(id, GL_COMPILE_STATUS, &status);
  if (status != GL_TRUE)
  {
    char buffer[512];
    if (type == GL_VERTEX_SHADER)
      std::cerr << "Vertex shader:" << std::endl;
    else if (type == GL_FRAGMENT_SHADER)
      std::cerr << "Fragment shader:" << std::endl;
    std::cerr << str << std::endl << std::endl;
    glGetShaderInfoLog(id, 512, NULL, buffer);
    std::cerr << "Error: " << std::endl << buffer << std::endl;
  }
  return id;
}
