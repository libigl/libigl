// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "texture_from_png.h"

#include "../opengl/report_gl_error.h"
#include <stb_image.h>

IGL_INLINE bool igl::png::texture_from_png(const std::string png_file, const bool flip, GLuint & id)
{
  int width,height,n;
  unsigned char *data = stbi_load(png_file.c_str(), &width, &height, &n, 4);
  if(data == NULL) {
    return false;
  }

  // Why do I need to flip?
  /*if(flip)
  {
    yimg.flip();
  }*/

  glGenTextures(1, &id);
  glBindTexture(GL_TEXTURE_2D, id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(
    GL_TEXTURE_2D, 0, GL_RGB,
    width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
  glBindTexture(GL_TEXTURE_2D, 0);

  stbi_image_free(data);

  return true;
}

IGL_INLINE bool igl::png::texture_from_png(const std::string png_file, GLuint & id)
{
  return texture_from_png(png_file,false,id);
}


