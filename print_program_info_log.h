#ifndef  IGL_PRINT_PROGRAM_INFO_LOG_H
#define IGL_PRINT_PROGRAM_INFO_LOG_H

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
namespace igl
{
  // Inputs:
  //   obj  OpenGL index of program to print info log about
  inline void print_program_info_log(const GLuint obj);
}

// Implmentation
#include <cstdio>
// Copyright Denis Kovacs 4/10/08
inline void igl::print_program_info_log(const GLuint obj)
{
  GLint infologLength = 0;
  GLint charsWritten  = 0;
  char *infoLog;
  
  glGetProgramiv(obj, GL_INFO_LOG_LENGTH,&infologLength);
  
  if (infologLength > 0)
  {
    infoLog = (char *)malloc(infologLength);
    glGetProgramInfoLog(obj, infologLength, &charsWritten, infoLog);
    printf("%s\n",infoLog);
    free(infoLog);
  }
}
#endif 
