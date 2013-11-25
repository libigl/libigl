// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "create_shader_program.h"
#ifndef IGL_NO_OPENGL

#include "load_shader.h"
#include "print_program_info_log.h"
#include <cstdio>

IGL_INLINE bool igl::create_shader_program(
  const std::string vert_source,
  const std::string frag_source,
  const std::map<std::string,GLuint> attrib,
  GLuint & id)
{
  if(vert_source == "" && frag_source == "")
  {
    fprintf(
      stderr,
      "Error: create_shader_program() could not create shader program,"
      " both .vert and .frag source given were empty\n");
    return false;
  }

  // create program
  id = glCreateProgram();
  if(id == 0)
  {
    fprintf(
      stderr,
      "Error: create_shader_program() could not create shader program.\n");
    return false;
  }

  if(vert_source != "")
  {
    // load vertex shader
    GLuint v = igl::load_shader(vert_source.c_str(),GL_VERTEX_SHADER);
    if(v == 0)
    {
      return false;
    }
    glAttachShader(id,v);
  }

  if(frag_source != "")
  {
    // load fragment shader
    GLuint f = igl::load_shader(frag_source.c_str(),GL_FRAGMENT_SHADER);
    if(f == 0)
    {
      return false;
    }
    glAttachShader(id,f);
  }

  // loop over attributes
  for(
    std::map<std::string,GLuint>::const_iterator ait = attrib.begin();
    ait != attrib.end();
    ait++)
  {
    glBindAttribLocation(
      id,
      (*ait).second,
      (*ait).first.c_str());
  }
  // Link program
  glLinkProgram(id);

  // print log if any
  igl::print_program_info_log(id);

  return true;
}
#endif
