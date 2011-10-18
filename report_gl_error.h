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
#include <string>

namespace igl
{
  // Print last OpenGL error to stderr prefixed by specified id string
  // Inputs:
  //   id   string to appear before any error msgs
  // Returns result of glGetError() 
  inline GLenum report_gl_error(const std::string id = std::string(""));
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

#endif
