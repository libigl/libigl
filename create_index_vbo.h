#ifndef IGL_CREATE_INDEX_VBO_H
#define IGL_CREATE_INDEX_VBO_H
// NOTE: It wouldn't be so hard to template this using Eigen's templates

#include <Eigen/Core>

#if __APPLE__
#  include <OpenGL/gl.h>
#else
#  ifdef _WIN32
#    define NOMINMAX
#    include <Windows.h>
#    undef NOMINMAX
#  endif
#  include <GL/gl.h>
#endif

// Create a VBO (Vertex Buffer Object) for a list of indices:
// GL_ELEMENT_ARRAY_BUFFER_ARB for the triangle indices (F)
namespace igl
{

  // Inputs:
  //   F  #F by 3 eigen Matrix of face (triangle) indices
  // Outputs:
  //   F_vbo_id  buffer id for face indices
  //
  inline void create_index_vbo(
    const Eigen::MatrixXi & F,
    GLuint & F_vbo_id);
}

// Implementation

// http://www.songho.ca/opengl/gl_vbo.html#create
inline void igl::create_index_vbo(
  const Eigen::MatrixXi & F,
  GLuint & F_vbo_id)
{
  // Generate Buffers
  glGenBuffersARB(1,&F_vbo_id);
  // Bind Buffers
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB,F_vbo_id);
  // Copy data to buffers
  // We expect a matrix with each vertex position on a row, we then want to
  // pass this data to OpenGL reading across rows (row-major)
  if(F.Options & Eigen::RowMajor)
  {
    glBufferDataARB(
      GL_ELEMENT_ARRAY_BUFFER_ARB,
      sizeof(int)*F.size(),
      F.data(),
      GL_STATIC_DRAW_ARB);
  }else
  {
    // Create temporary copy of transpose
    Eigen::MatrixXi FT = F.transpose();
    // If its column major then we need to temporarily store a transpose
    glBufferDataARB(
      GL_ELEMENT_ARRAY_BUFFER_ARB,
      sizeof(int)*F.size(),
      FT.data(),
      GL_STATIC_DRAW);
  }
  // bind with 0, so, switch back to normal pointer operation
  glBindBufferARB(GL_ELEMENT_ARRAY_BUFFER_ARB, 0);
}

#endif

