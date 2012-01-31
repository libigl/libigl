#include "per_face_normals.h"

IGL_INLINE void igl::per_face_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & N)
{
  N.resize(F.rows(),3);
  // loop over faces
  for(int i = 0; i < F.rows();i++)
  {
    float v1[3];
    v1[0] = V(F(i,1),0) - V(F(i,0),0);
    v1[1] = V(F(i,1),1) - V(F(i,0),1);
    v1[2] = V(F(i,1),2) - V(F(i,0),2);
    float v2[3];
    v2[0] = V(F(i,2),0) - V(F(i,0),0);
    v2[1] = V(F(i,2),1) - V(F(i,0),1);
    v2[2] = V(F(i,2),2) - V(F(i,0),2);
    N(i,0) = v1[1]*v2[2] - v1[2]*v2[1];
    N(i,1) = -(v1[0]*v2[2] - v1[2]*v2[0]);
    N(i,2) = v1[0]*v2[1] - v1[1]*v2[0];
    float length = sqrt(
      N(i,0)*N(i,0) +
      N(i,1)*N(i,1) +
      N(i,2)*N(i,2));
    N(i,0) /= length;
    N(i,1) /=  length;
    N(i,2) /=  length;
  }
}
