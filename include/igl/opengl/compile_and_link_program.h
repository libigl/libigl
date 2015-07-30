// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_COMPILE_AND_LINK_PROGRAM_H
#define IGL_OPENGL_COMPILE_AND_LINK_PROGRAM_H
#include "../igl_inline.h"
#include "OpenGL_convenience.h"
namespace igl
{
  namespace opengl
  {
    // Compile and link very simple vertex/fragment shaders
    //
    // Inputs:
    //   v_str  string of vertex shader contents
    //   f_str  string of fragment shader contents
    // Returns id of program
    //
    // Known bugs: this seems to duplicate `create_shader_program` with less
    // functionality.
    IGL_INLINE GLuint compile_and_link_program(
      const char * v_str, const char * f_str);
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "compile_and_link_program.cpp"
#endif
#endif
