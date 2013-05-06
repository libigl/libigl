#ifndef IGL_CREATE_SHADER_PROGRAM_H
#define IGL_CREATE_SHADER_PROGRAM_H
#include "igl_inline.h"
#include <string>
#include <map>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#elif defined(_WIN32)
#  define NOMINMAX
#  include <Windows.h>
#  undef NOMINMAX
#  include <GL/glew.h>
#  include <GL/gl.h>
#else
#  define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
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
  IGL_INLINE bool create_shader_program(
    const std::string vert_source,
    const std::string frag_source,
    const std::map<std::string,GLuint> attrib,
    GLuint & id);
}

#ifdef IGL_HEADER_ONLY
#  include "create_shader_program.cpp"
#endif

#endif
