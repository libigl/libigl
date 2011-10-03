#ifndef IGL_CREATE_SHADER_PROGRAM_H
#define IGL_CREATE_SHADER_PROGRAM_H
#include <string>
#include <map>

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#else
#   include <GL/gl.h>
#endif

namespace igl
{
  // Create a shader program with a vertex and fragments shader loading from
  // source strings and vertex attributes assigned from a map before linking the
  // shaders to the program, making it ready to use with glUseProgram(id)
  // Inputs:
  //   vert_source  string containing source code of vertex shader
  //   frag_source  string containing source code of fragment shader
  //   attrib  map containing table of vertex attribute strings add their
  //   correspondingly ids (generated previously using glBindAttribLocation)
  // Outputs:
  //   id  index id of created shader, set to 0 on error
  // Returns true on success, false on error
  //
  // Note: Caller is responsible for making sure that current value of id is not
  // leaking a shader (since it will be overwritten)
  //
  // See also: destroy_shader_program
  inline bool create_shader_program(
    const std::string vert_source,
    const std::string frag_source,
    const std::map<std::string,GLuint> attrib,
    GLuint & id);
}

// Implementation
#include "load_shader.h"
#include "print_program_info_log.h"
#include <cstdio>

inline bool igl::create_shader_program(
  const std::string vert_source,
  const std::string frag_source,
  const std::map<std::string,GLuint> attrib,
  GLuint & id)
{
  if(vert_source == "" && frag_source == "")
  {
    fprintf(
      stderr,
      "Error: create_shader_program() could not create shader program,"
      " both .vert and .frag source given were empty\n");
    return false;
  }

  // create program
  id = glCreateProgram();
  if(id == 0)
  {
    fprintf(
      stderr,
      "Error: create_shader_program() could not create shader program.\n");
    return false;
  }

  if(vert_source != "")
  {
    // load vertex shader
    GLuint v = igl::load_shader(vert_source.c_str(),GL_VERTEX_SHADER);
    if(v == 0)
    {
      return false;
    }
    glAttachShader(id,v);
  }

  if(frag_source != "")
  {
    // load fragment shader
    GLuint f = igl::load_shader(frag_source.c_str(),GL_FRAGMENT_SHADER);
    if(f == 0)
    {
      return false;
    }
    glAttachShader(id,f);
  }

  // loop over attributes
  for(
    std::map<std::string,GLuint>::const_iterator ait = attrib.begin();
    ait != attrib.end();
    ait++)
  {
    glBindAttribLocation(
      id,
      (*ait).second,
      (*ait).first.c_str());
  }
  // Link program
  glLinkProgram(id);
  // print log if any
  igl::print_program_info_log(id);

  return true;
}
#endif
