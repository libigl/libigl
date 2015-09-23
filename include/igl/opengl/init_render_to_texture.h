// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_INIT_RENDER_TO_TEXTURE_H
#define IGL_OPENGL_INIT_RENDER_TO_TEXTURE_H
#include "../igl_inline.h"
#include "OpenGL_convenience.h"
#include <cstdlib>
namespace igl
{
  namespace opengl
  {
    // Create a texture+framebuffer+depthbuffer triplet bound for rendering into
    // the texture;
    //
    // Inputs:
    //   width  image width
    //   height image height
    // Outputs:
    //   tex_id  id of the texture
    //   fbo_id  id of the frame buffer object
    //   dfbo_id  id of the depth frame buffer object
    IGL_INLINE void init_render_to_texture(
      const size_t width,
      const size_t height,
      GLuint & tex_id,
      GLuint & fbo_id,
      GLuint & dfbo_id);
  }
}
#ifndef IGL_STATIC_LIBRARY
#  include "init_render_to_texture.cpp"
#endif
#endif
