// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "texture_from_png.h"

#include "../opengl/report_gl_error.h"
#include <YImage.hpp>

IGL_INLINE bool igl::png::texture_from_png(const std::string png_file, GLuint & id)
{
  YImage yimg;
  if(!yimg.load(png_file.c_str()))
  {
    return false;
  }
  // Why do I need to flip?
  //yimg.flip();
  glGenTextures(1, &id);
  glBindTexture(GL_TEXTURE_2D, id);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexImage2D(
    GL_TEXTURE_2D, 0, GL_RGB,
    yimg.width(), yimg.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, yimg.data());
  glBindTexture(GL_TEXTURE_2D, 0);
  return true;
}


IGL_INLINE bool igl::png::texture_from_png(
  const std::string png_file,
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& R,
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& G,
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& B,
  Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& A
)
{
  YImage yimg;
  if(!yimg.load(png_file.c_str()))
  {
    return false;
  }

  R.resize(yimg.height(),yimg.width());
  G.resize(yimg.height(),yimg.width());
  B.resize(yimg.height(),yimg.width());
  A.resize(yimg.height(),yimg.width());

  for (unsigned j=0; j<yimg.height(); ++j)
  {
    for (unsigned i=0; i<yimg.width(); ++i)
    {
      R(i,j) = yimg.at(yimg.width()-1-i,yimg.height()-1-j).r;
      G(i,j) = yimg.at(yimg.width()-1-i,yimg.height()-1-j).g;
      B(i,j) = yimg.at(yimg.width()-1-i,yimg.height()-1-j).b;
      //1A(i,j) = yimg.at(yimg.width()-1-i,yimg.height()-1-j).a;
    }
  }

  return true;
}
