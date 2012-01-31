#include "per_vertex_normals.h"

#include "per_face_normals.h"
#include "normalize_rows.h"

IGL_INLINE void igl::per_vertex_normals(
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
