#include "create_vector_vbo.h"

#include <cassert>

// http://www.songho.ca/opengl/gl_vbo.html#create
template <typename T>
IGL_INLINE void igl::create_vector_vbo(
  const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & V,
  GLuint & V_vbo_id)
{
  //// Expects that input is list of 3D vectors along rows
  //assert(V.cols() == 3);

  // Generate Buffers
  glGenBuffers(1,&V_vbo_id);
  // Bind Buffers
  glBindBuffer(GL_ARRAY_BUFFER,V_vbo_id);
  // Copy data to buffers
  // We expect a matrix with each vertex position on a row, we then want to
  // pass this data to OpenGL reading across rows (row-major)
  if(V.Options & Eigen::RowMajor)
  {
    glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(T)*V.size(),
      V.data(),
      GL_STATIC_DRAW);
  }else
  {
    // Create temporary copy of transpose
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> VT = V.transpose();
    // If its column major then we need to temporarily store a transpose
    glBufferData(
      GL_ARRAY_BUFFER,
      sizeof(T)*V.size(),
      VT.data(),
      GL_STATIC_DRAW);
  }
  // bind with 0, so, switch back to normal pointer operation
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
