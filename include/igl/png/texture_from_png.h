// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_TEXTURE_FROM_PNG_H
#define IGL_TEXTURE_FROM_PNG_H
#ifndef IGL_NO_OPENGL
#include "../igl_inline.h"

#include "../OpenGL_convenience.h"

#include <string>

namespace igl
{
  // Read an image from a .png file and use it as a texture
  //
  // Input:
  //  png_file  path to .png file
  // Output:
  //  id  of generated openGL texture
  // Returns true on success, false on failure
  IGL_INLINE bool texture_from_png(const std::string png_file, GLuint & id);
}

#ifndef IGL_STATIC_LIBRARY
#  include "texture_from_png.cpp"
#endif

#endif
#endif

