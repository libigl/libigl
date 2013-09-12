#include "report_gl_error.h"
#ifndef IGL_NO_OPENGL

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
