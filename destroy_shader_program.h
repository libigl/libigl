#ifndef IGL_DESTROY_SHADER_PROGRAM_H
#define IGL_DESTROY_SHADER_PROGRAM_H

#ifdef __APPLE__
#   include <OpenGL/gl.h>
#else
#   include <GL/gl.h>
#endif

namespace igl
{
  // Properly destroy a shader program. Detach and delete each of its shaders
  // and delete it
  // Inputs:
  //   id  index id of created shader, set to 0 on error
  // Returns true on success, false on error
  // 
  // Note: caller is responsible for making sure he doesn't foolishly continue
  // to use id as if it still contains a program
  // 
  // See also: create_shader_program
  inline bool destroy_shader_program(const GLuint id);
}

// Implementation
inline bool igl::destroy_shader_program(const GLuint id)
{
  // Don't try to destroy id == 0 (no shader program)
  if(id == 0)
  {
    fprintf(stderr,"Error: destroy_shader_program() id = %d"
      " but must should be positive\n",id);
    return false;
  }
  // Get each attached shader one by one and detach and delete it
  GLsizei count;
  // shader id
  GLuint s;
  do
  {
    // Try to get at most *1* attached shader
    glGetAttachedShaders(id,1,&count,&s);
    // Check that we actually got *1*
    if(count == 1)
    {
      // Detach and delete this shader
      glDetachShader(id,s);
      glDeleteShader(s);
    }
  }while(count > 0);
  // Now that all of the shaders are gone we can just delete the program
  glDeleteProgram(id);
  return true;
}

#endif
