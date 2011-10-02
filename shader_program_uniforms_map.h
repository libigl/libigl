#ifndef IGL_SHADER_PROGRAM_UNIFORMS_MAP_H
#define IGL_SHADER_PROGRAM_UNIFORMS_MAP_H
#include <string>
#include <map>

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#else
#   include <GL/gl.h>
#endif

namespace igl
{
  // Builds a map of *active* uniform names as strings to their respective
  // indices (NOT locations) as GLuint.
  // Input:
  //   id  index id of the program to query 
  // Output:
  //   uniforms  map of names to indices
  // Returns true on success, false on errors
  void shader_program_uniforms_map(
    const GLuint id, 
    std::map<std::string,GLint> & uniforms);
}

// Implementation
#include "verbose.h"
#include "report_gl_error.h"
#include "uniform_type_to_string.h"

void igl::shader_program_uniforms_map(
  const GLuint id, 
  std::map<std::string,GLint> & uniforms)
{
  // empty the map of previous contents
  uniforms.clear();

  // get number of active uniforms
  GLint n = 200;
  glGetProgramiv(id,GL_ACTIVE_UNIFORMS,&n);

  // get max uniform name length
  GLint max_name_length;
  glGetProgramiv(id,GL_ACTIVE_UNIFORM_MAX_LENGTH,&max_name_length);

  // buffer for name
  GLchar * name = new GLchar[max_name_length];
  // buffer for length
  GLsizei length = 100;
  GLenum type;
  GLint size;

  // loop over active uniforms getting each's name
  for(GLuint u = 0;u < n;u++)
  {
    // I have no idea why glGetActiveUniformName doesn't work but
    // glGetActiveUniform does...
    //glGetActiveUniformName(id,u,max_name_length,&length,name);
    glGetActiveUniform(id,u,max_name_length,&length,&size,&type,name);
    // insert into map
    uniforms[string(name)] = u;
    verbose("%s --> index: %d size: %d type: %s\n",
      name,
      (int)u,
      size,
      uniform_type_to_string(type).c_str());
  }

  delete[] name;
}
#endif
