#ifndef IGL_REPORT_GL_ERROR
#define IGL_REPORT_GL_ERROR

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#   include <OpenGL/glu.h>
#else
#   include <GL/gl.h>
#   include <GL/glu.h>
#endif

#include <cstdio>

namespace igl
{
  // Print last OpenGL error to stderr
  // Returns result of glGetError()
  inline GLenum report_gl_error();
}

// Implementation
#include "verbose.h"

inline GLenum igl::report_gl_error()
{
  GLenum err = glGetError();
  if(GL_NO_ERROR != err)
  {
    verbose("GL_ERROR: ");
    fprintf(stderr,"%s\n",gluErrorString(err));
  }
  return err;
}

#endif
