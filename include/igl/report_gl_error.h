#ifndef IGL_REPORT_GL_ERROR_H
#define IGL_REPORT_GL_ERROR_H
#ifndef IGL_NO_OPENGL
#include "igl_inline.h"

#include "OpenGL_convenience.h"

#include <cstdio>
#include <string>

namespace igl
{
  // Print last OpenGL error to stderr prefixed by specified id string
  // Inputs:
  //   id   string to appear before any error msgs
  // Returns result of glGetError() 
  IGL_INLINE GLenum report_gl_error(const std::string id);
  // No prefix
  IGL_INLINE GLenum report_gl_error();
}

#ifdef IGL_HEADER_ONLY
#  include "report_gl_error.cpp"
#endif

#endif
#endif
