// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "OpenGL_shader.h"

#ifndef __APPLE__
#  define GLEW_STATIC
#  include <GL/glew.h>
#endif

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   ifdef _WIN32
#       include <windows.h>
#   endif
#   include <GL/gl.h>
#endif

#include <iostream>
#include <fstream>

IGL_INLINE bool igl::viewer::OpenGL_shader::init_from_files(
  const std::string &vertex_shader_filename,
  const std::string &fragment_shader_filename,
  const std::string &fragment_data_name,
  const std::string &geometry_shader_filename,
  int geometry_shader_max_vertices)
{
  auto file_to_string = [](const std::string &filename)->std::string
  {
    std::ifstream t(filename);
    return std::string((std::istreambuf_iterator<char>(t)),
                        std::istreambuf_iterator<char>());
  };

  return init(
    file_to_string(vertex_shader_filename),
    file_to_string(fragment_shader_filename),
    fragment_data_name,
    file_to_string(geometry_shader_filename),
    geometry_shader_max_vertices
 );
}

IGL_INLINE bool igl::viewer::OpenGL_shader::init(
  const std::string &vertex_shader_string,
  const std::string &fragment_shader_string,
  const std::string &fragment_data_name,
  const std::string &geometry_shader_string,
  int geometry_shader_max_vertices)
{
  using namespace std;
  vertex_shader = create_shader_helper(GL_VERTEX_SHADER, vertex_shader_string);
  geometry_shader = create_shader_helper(GL_GEOMETRY_SHADER, geometry_shader_string);
  fragment_shader = create_shader_helper(GL_FRAGMENT_SHADER, fragment_shader_string);

  if (!vertex_shader || !fragment_shader)
    return false;

  program_shader = glCreateProgram();

  glAttachShader(program_shader, vertex_shader);
  glAttachShader(program_shader, fragment_shader);

  if (geometry_shader)
  {
    glAttachShader(program_shader, geometry_shader);

    /* This covers only basic cases and may need to be modified */
    glProgramParameteri(program_shader, GL_GEOMETRY_INPUT_TYPE, GL_TRIANGLES);
    glProgramParameteri(program_shader, GL_GEOMETRY_OUTPUT_TYPE, GL_TRIANGLES);
    glProgramParameteri(program_shader, GL_GEOMETRY_VERTICES_OUT, geometry_shader_max_vertices);
  }

  glBindFragDataLocation(program_shader, 0, fragment_data_name.c_str());
  glLinkProgram(program_shader);

  GLint status;
  glGetProgramiv(program_shader, GL_LINK_STATUS, &status);

  if (status != GL_TRUE)
  {
    char buffer[512];
    glGetProgramInfoLog(program_shader, 512, NULL, buffer);
    cerr << "Linker error: " << endl << buffer << endl;
    program_shader = 0;
    return false;
  }

  return true;
}

IGL_INLINE void igl::viewer::OpenGL_shader::bind()
{
  glUseProgram(program_shader);
}

IGL_INLINE GLint igl::viewer::OpenGL_shader::attrib(const std::string &name) const
{
  return glGetAttribLocation(program_shader, name.c_str());
}

IGL_INLINE GLint igl::viewer::OpenGL_shader::uniform(const std::string &name) const
{
  return glGetUniformLocation(program_shader, name.c_str());
}

IGL_INLINE GLint igl::viewer::OpenGL_shader::bindVertexAttribArray(
  const std::string &name, GLuint bufferID, const Eigen::MatrixXf &M, bool refresh) const
{
  GLint id = attrib(name);
  if (id < 0)
    return id;
  if (M.size() == 0)
  {
    glDisableVertexAttribArray(id);
    return id;
  }
  glBindBuffer(GL_ARRAY_BUFFER, bufferID);
  if (refresh)
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*M.size(), M.data(), GL_DYNAMIC_DRAW);
  glVertexAttribPointer(id, M.rows(), GL_FLOAT, GL_FALSE, 0, 0);
  glEnableVertexAttribArray(id);
  return id;
}

IGL_INLINE void igl::viewer::OpenGL_shader::free()
{
  if (program_shader)
  {
    glDeleteProgram(program_shader);
    program_shader = 0;
  }
  if (vertex_shader)
  {
    glDeleteShader(vertex_shader);
    vertex_shader = 0;
  }
  if (fragment_shader)
  {
    glDeleteShader(fragment_shader);
    fragment_shader = 0;
  }
  if (geometry_shader)
  {
    glDeleteShader(geometry_shader);
    geometry_shader = 0;
  }
}

IGL_INLINE GLuint igl::viewer::OpenGL_shader::create_shader_helper(GLint type, const std::string &shader_string)
{
  using namespace std;
  if (shader_string.empty())
    return (GLuint) 0;

  GLuint id = glCreateShader(type);
  const char *shader_string_const = shader_string.c_str();
  glShaderSource(id, 1, &shader_string_const, NULL);
  glCompileShader(id);

  GLint status;
  glGetShaderiv(id, GL_COMPILE_STATUS, &status);

  if (status != GL_TRUE)
  {
    char buffer[512];
    if (type == GL_VERTEX_SHADER)
      cerr << "Vertex shader:" << endl;
    else if (type == GL_FRAGMENT_SHADER)
      cerr << "Fragment shader:" << endl;
    else if (type == GL_GEOMETRY_SHADER)
      cerr << "Geometry shader:" << endl;
    cerr << shader_string << endl << endl;
    glGetShaderInfoLog(id, 512, NULL, buffer);
    cerr << "Error: " << endl << buffer << endl;
    return (GLuint) 0;
  }

  return id;
}
