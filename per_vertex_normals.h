#ifndef IGL_PER_VERTEX_NORMALS_H
#define IGL_PER_VERTEX_NORMALS_H
#include <Eigen/Core>
// Note: So for this only computes normals per vertex as uniformly weighted
// averages of incident triangle normals. It would be nice to support more or
// all of the methods here:
// "A comparison of algorithms for vertex normal computation"
namespace igl
{
  // Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  // Output:
  //   N  #V by 3 eigen Matrix of mesh vertex 3D normals
  void per_vertex_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    Eigen::MatrixXd & N);
}

// Implementation
#include "per_face_normals.h"
#include "normalize_rows.h"

void igl::per_vertex_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  Eigen::MatrixXd PFN;
  igl::per_face_normals(V,F,PFN);

  // Resize for output
  N.resize(V.rows(),3);
  // loop over vertices, setting normalize to 0
  for(int i = 0; i < N.rows();i++)
  {
    N(i,0) = 0;
    N(i,1) = 0;
    N(i,2) = 0;
  }

  // loop over faces
  for(int i = 0; i < F.rows();i++)
  {
    // throw normal at each corner
    for(int j = 0; j < 3;j++)
    {
      N(F(i,j),0) += PFN(i,0);
      N(F(i,j),1) += PFN(i,1);
      N(F(i,j),2) += PFN(i,2);
    }
  }
  // normalize each row
  igl::normalize_rows(N,N);
}

#endif
