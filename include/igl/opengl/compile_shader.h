// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_COMPILE_SHADER_H
#define IGL_OPENGL_COMPILE_SHADER_H
#include "OpenGL_convenience.h"
#include "../igl_inline.h"
namespace igl
{
  namespace opengl
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
}
#ifndef IGL_STATIC_LIBRARY
#  include "compile_shader.cpp"
#endif
#endif
