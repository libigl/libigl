#ifndef IGL_PER_CORNER_NORMALS_H
#define IGL_PER_CORNER_NORMALS_H
#include <Eigen/Core>
#endif

namespace igl
{
  // Compute vertex normals via vertex position list, face list
  // Inputs:
  //   V  #V by 3 eigen Matrix of mesh vertex 3D positions
  //   F  #F by 3 eigne Matrix of face (triangle) indices
  //   corner_threshold  threshold in degrees on sharp angles
  // Output:
  //   CN  #F*3 by 3 eigen Matrix of mesh vertex 3D normals, where the normal
  //     for corner F(i,j) is at CN(i*3+j,:) 
  inline void per_corner_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const double corner_threshold,
    Eigen::MatrixXd & CN);
  // Other Inputs:
  //   FN  #F by 3 eigen Matrix of face normals
  //   VF  map from vertices to list of incident faces
  inline void per_corner_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const std::vector<std::vector<int> >& VF,
    const double corner_threshold,
    Eigen::MatrixXd & CN);
}

// Implementation
#include "vf.h"
#include "per_face_normals.h"
#include "PI.h"

inline void igl::per_corner_normals(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const double corner_threshold,
  Eigen::MatrixXd & CN)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;
  MatrixXd FN;
  per_face_normals(V,F,FN);
  vector<vector<int> > VF,VFi;
  vf(V,F,VF,VFi);
  return per_corner_normals(V,F,FN,VF,corner_threshold,CN);
}

inline void igl::per_corner_normals(
    const Eigen::MatrixXd & V,
    const Eigen::MatrixXi & F,
    const Eigen::MatrixXd & FN,
    const std::vector<std::vector<int> >& VF,
    const double corner_threshold,
    Eigen::MatrixXd & CN)
{
  using namespace igl;
  using namespace Eigen;
  using namespace std;

  // number of faces
  const int m = F.rows();
  // initialize output to zero
  CN = MatrixXd::Constant(m*3,3,0);

  // loop over faces
  for(size_t i = 0;i<m;i++)
  {
    // Normal of this face
    Vector3d fn = FN.row(i);
    // loop over corners
    for(size_t j = 0;j<3;j++)
    {
      std::vector<int> incident_faces = VF[F(i,j)];
      // loop over faces sharing vertex of this corner
      for(int k = 0;k<(int)incident_faces.size();k++)
      {
        Vector3d ifn = FN.row(incident_faces[k]);
        // dot product between face's normal and other face's normal
        double dp = fn.dot(ifn);
        // if difference in normal is slight then add to average
        if(dp > cos(corner_threshold*PI/180))
        {
          // add to running sum
          CN.row(i*3+j) += ifn;
        // else ignore
        }else
        {
        }
      }
      // normalize to take average
      double length = CN.row(i*3+j).norm();
      CN.row(i*3+j) /= length;
    }
  }
}

