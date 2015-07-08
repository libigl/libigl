// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Wenzel Jacob <wenzel@inf.ethz.ch>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_VIEWER_OPENGL_SHADER_H
#define IGL_VIEWER_OPENGL_SHADER_H

#include <igl/igl_inline.h>
#include <string>
#include <Eigen/Core>

#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#  undef DrawText
#endif

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

namespace igl
{
namespace viewer
{

// This class wraps an OpenGL program composed of three shaders
// TODO: write documentation

class OpenGL_shader
{
public:
  typedef unsigned int GLuint;
  typedef int GLint;

  GLuint vertex_shader;
  GLuint fragment_shader;
  GLuint geometry_shader;
  GLuint program_shader;

  IGL_INLINE OpenGL_shader() : vertex_shader(0), fragment_shader(0),
    geometry_shader(0), program_shader(0) { }

  // Create a new shader from the specified source strings
  IGL_INLINE bool init(const std::string &vertex_shader_string,
    const std::string &fragment_shader_string,
    const std::string &fragment_data_name,
    const std::string &geometry_shader_string = "",
    int geometry_shader_max_vertices = 3);

  // Create a new shader from the specified files on disk
  IGL_INLINE bool init_from_files(const std::string &vertex_shader_filename,
    const std::string &fragment_shader_filename,
    const std::string &fragment_data_name,
    const std::string &geometry_shader_filename = "",
    int geometry_shader_max_vertices = 3);

  // Select this shader for subsequent draw calls
  IGL_INLINE void bind();

  // Release all OpenGL objects
  IGL_INLINE void free();

  // Return the OpenGL handle of a named shader attribute (-1 if it does not exist)
  IGL_INLINE GLint attrib(const std::string &name) const;

  // Return the OpenGL handle of a uniform attribute (-1 if it does not exist)
  IGL_INLINE GLint uniform(const std::string &name) const;

  // Bind a per-vertex array attribute and refresh its contents from an Eigen amtrix
  IGL_INLINE GLint bindVertexAttribArray(const std::string &name, GLuint bufferID,
    const Eigen::MatrixXf &M, bool refresh) const;

  IGL_INLINE GLuint create_shader_helper(GLint type, const std::string &shader_string);

};

}
}

#ifndef IGL_STATIC_LIBRARY
#  include "OpenGL_shader.cpp"
#endif

#endif
