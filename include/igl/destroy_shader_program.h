// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_DESTROY_SHADER_PROGRAM_H
#define IGL_DESTROY_SHADER_PROGRAM_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

namespace igl
{
  // Properly destroy a shader program. Detach and delete each of its shaders
  // and delete it
  // Inputs:
  //   id  index id of created shader, set to 0 on error
  // Returns true on success, false on error
  // 
  // Note: caller is responsible for making sure he doesn't foolishly continue
  // to use id as if it still contains a program
  // 
  // See also: create_shader_program
  IGL_INLINE bool destroy_shader_program(const GLuint id);
}

#ifndef IGL_STATIC_LIBRARY
#  include "destroy_shader_program.cpp"
#endif

#endif
#endif
