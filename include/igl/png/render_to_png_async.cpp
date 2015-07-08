// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "render_to_png_async.h"
#include <YImage.hpp>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#endif

  
static IGL_INLINE bool render_to_png_async_helper(
  YImage * img,
  const std::string png_file,
  const bool alpha,
  const bool fast)
{

  img->flip();
  const int width = img->width();
  const int height = img->height();
  if(!alpha)
  {
    for(int i = 0;i<width;i++)
    for(int j = 0;j<height;j++)
    {
      img->at(i,j).a = 255;
    }
  }

  bool ret = img->save(png_file.c_str(),fast);
  delete img;
  return ret;
}

IGL_INLINE std::thread igl::png::render_to_png_async(
  const std::string png_file,
  const int width,
  const int height,
  const bool alpha,
  const bool fast)
{
  // Part that should serial
  YImage * img = new YImage();
  img->resize(width,height);
  glReadPixels(
    0,
    0,
    width,
    height,
    GL_RGBA,
    GL_UNSIGNED_BYTE,
    img->data());
  // Part that should be asynchronous  
  std::thread t(render_to_png_async_helper,img,png_file,alpha,fast);
  t.detach();
  return t;
}
