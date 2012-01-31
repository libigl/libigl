#ifndef IGL_REPORT_GL_ERROR_H
#define IGL_REPORT_GL_ERROR_H

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#  include <OpenGL/glu.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#  include <GL/glu.h>
#endif

#include <cstdio>
#include <string>

namespace igl
{
  // Print last OpenGL error to stderr prefixed by specified id string
  // Inputs:
  //   id   string to appear before any error msgs
  // Returns result of glGetError() 
  inline GLenum report_gl_error(const std::string id);
  // No prefix
  inline GLenum report_gl_error();
}

// Implementation
#include "verbose.h"

inline GLenum igl::report_gl_error(const std::string id)
{
  GLenum err = glGetError();
  if(GL_NO_ERROR != err)
  {
    verbose("GL_ERROR: ");
    fprintf(stderr,"%s%s\n",id.c_str(),gluErrorString(err));
  }
  return err;
}

inline GLenum igl::report_gl_error()
{
  return igl::report_gl_error(std::string(""));
}

#endif
