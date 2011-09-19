#ifndef IGL_CREATE_VECTOR_VBO
#define IGL_CREATE_VECTOR_VBO
// NOTE: It wouldn't be so hard to template this using Eigen's templates

#include <Eigen/Core>

#if __APPLE__
#  include <OpenGL/gl.h>
#else
#  include <GL/gl.h>
#endif

// Create a VBO (Vertex Buffer Object) for a list of vectors:
// GL_ARRAY_BUFFER_ARB for the vectors (V)
namespace igl
{

  // Inputs:
  //   V  #V by 3 eigen Matrix of vertex 3D positions
  // Outputs:
  //   V_vbo_id  buffer id for vectors
  //
  void create_vector_vbo(
    const Eigen::MatrixXd & V,
    GLuint & V_vbo_id);
}

// Implementation
#include <cassert>

// http://www.songho.ca/opengl/gl_vbo.html#create
void igl::create_vector_vbo(
  const Eigen::MatrixXd & V,
  GLuint & V_vbo_id)
{
  // Expects that input is list of 3D vectors along rows
  assert(V.cols() == 3);

  // Generate Buffers
  glGenBuffersARB(1,&V_vbo_id);
  // Bind Buffers
  glBindBufferARB(GL_ARRAY_BUFFER_ARB,V_vbo_id);
  // Copy data to buffers
  // We expect a matrix with each vertex position on a row, we then want to
  // pass this data to OpenGL reading across rows (row-major)
  if(V.Options & Eigen::RowMajor)
  {
    glBufferDataARB(
      GL_ARRAY_BUFFER_ARB,
      sizeof(double)*V.size(),
      V.data(),
      GL_STATIC_DRAW_ARB);
  }else
  {
    // Create temporary copy of transpose
    Eigen::MatrixXd VT = V.transpose();
    // If its column major then we need to temporarily store a transpose
    glBufferDataARB(
      GL_ARRAY_BUFFER_ARB,
      sizeof(double)*V.size(),
      VT.data(),
      GL_STATIC_DRAW_ARB);
  }
  // bind with 0, so, switch back to normal pointer operation
  glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);
}

#endif

