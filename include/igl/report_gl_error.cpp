// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "report_gl_error.h"
#ifndef IGL_NO_OPENGL

#include <cstdio>
#include "verbose.h"

IGL_INLINE GLenum igl::report_gl_error(const std::string id)
{
  GLenum err = glGetError();
  if(GL_NO_ERROR != err)
  {
    verbose("GL_ERROR: ");
    fprintf(stderr,"%s%s\n",id.c_str(),gluErrorString(err));
  }
  return err;
}

IGL_INLINE GLenum igl::report_gl_error()
{
  return igl::report_gl_error(std::string(""));
}
#endif
