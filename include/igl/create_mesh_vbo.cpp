#include "create_mesh_vbo.h"

#include "create_vector_vbo.h"
#include "create_index_vbo.h"

// http://www.songho.ca/opengl/gl_vbo.html#create
IGL_INLINE void igl::create_mesh_vbo(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  GLuint & V_vbo_id,
  GLuint & F_vbo_id)
{
  // Create VBO for vertex position vectors
  create_vector_vbo(V,V_vbo_id);
  // Create VBO for face index lists
  create_index_vbo(F,F_vbo_id);
}

// http://www.songho.ca/opengl/gl_vbo.html#create
IGL_INLINE void igl::create_mesh_vbo(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::MatrixXd & N,
  GLuint & V_vbo_id,
  GLuint & F_vbo_id,
  GLuint & N_vbo_id)
{
  // Create VBOs for faces and vertices
  create_mesh_vbo(V,F,V_vbo_id,F_vbo_id);
  // Create VBO for normal vectors
  create_vector_vbo(N,N_vbo_id);
}

#ifndef IGL_HEADER_ONLY
// Explicit template specialization
#endif
