#ifndef IGL_CREATE_VECTOR_VBO_H
#define IGL_CREATE_VECTOR_VBO_H
#include "igl_inline.h"
// NOTE: It wouldn't be so hard to template this using Eigen's templates

#include <Eigen/Core>

#ifdef __APPLE__
#  include <OpenGL/gl.h>
#elif defined(_WIN32)
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#else
#define GL_GLEXT_PROTOTYPES
#  include <GL/gl.h>
#  include <GL/glext.h>
#endif

// Create a VBO (Vertex Buffer Object) for a list of vectors:
// GL_ARRAY_BUFFER for the vectors (V)
namespace igl
{

  // Templates:
  //   T  should be a eigen matrix primitive type like int or double
  // Inputs:
  //   V  m by n eigen Matrix of type T values
  // Outputs:
  //   V_vbo_id  buffer id for vectors
  //
  template <typename T>
  IGL_INLINE void create_vector_vbo(
    const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & V,
    GLuint & V_vbo_id);
}

#ifdef IGL_HEADER_ONLY
#  include "create_vector_vbo.cpp"
#endif

#endif
