#ifndef IGL_REPORT_GL_ERROR_H
#define IGL_REPORT_GL_ERROR_H
#include "igl_inline.h"

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
  IGL_INLINE GLenum report_gl_error(const std::string id);
  // No prefix
  IGL_INLINE GLenum report_gl_error();
}

#ifdef IGL_HEADER_ONLY
#  include "report_gl_error.cpp"
#endif

#endif
