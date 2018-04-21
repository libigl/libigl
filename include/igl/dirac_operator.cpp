#include <igl/igl_inline.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <igl/doublearea.h>
#include <iostream>
#include "dirac_operator.h"
template <typename DerivedV, typename DerivedF, typename Scalar>
IGL_INLINE void igl::dirac_operator(
const Eigen::MatrixBase<DerivedV> & V, 
const Eigen::MatrixBase<DerivedF> & F, 
Eigen::SparseMatrix<Scalar>& D,
Eigen::SparseMatrix<Scalar>& DA)
{

  auto quanternion3_matrix = [](Eigen::Matrix<Scalar, 1, 3> e) 
  {
    auto a = 0;
    auto b = e[0];
    auto c = e[1];
    auto d = e[2];
    Eigen::Matrix<Scalar, 4, 4> mat;
    mat << a, -b, -c, -d,
            b, a, -d, c,
            c, d, a, -b,
            d, -c, b, a;
    return mat;
  };

  Eigen::Matrix<Scalar, -1, 1> dblAf;
  igl::doublearea(V, F, dblAf);
  
  Eigen::Matrix<Scalar, -1, 1> dblAv(V.rows());
  dblAv.setZero();
  std::vector<Eigen::Triplet<Scalar>> IJV_D;
  std::vector<Eigen::Triplet<Scalar>> IJV_DA;
  
  for(int f=0; f<F.rows(); f++)
    for (int j=0; j<3; j++)
    {
      dblAv[F(f,j)] += dblAf[f] / 3.0;
    }
  
  IJV_D.reserve(4*4*3*F.rows());
  IJV_DA.reserve(4*4*3*F.rows());
  for(int i=0; i<F.rows(); i++)
  {
      for(int ind=0; ind<3; ind++)
      {
          int j = F(i, ind);
          int ind1 = F(i, (ind+1)%3);
          int ind2 = F(i, (ind+2)%3);
  
          auto e1 = V.row(ind1).template cast<Scalar>();
          auto e2 = V.row(ind2).template cast<Scalar>();
  
          auto mat = - quanternion3_matrix(e1 - e2);
  
          for(int ii=0; ii<4; ii++) {
              for (int jj=0; jj<4; jj++) {
                  IJV_D.push_back(Eigen::Triplet<Scalar>(4*i+ii, 4*j+jj, mat(ii,jj)/dblAf[i]));
                  IJV_DA.push_back(Eigen::Triplet<Scalar>(4*j+jj, 4*i+ii, mat(ii,jj)/dblAv[j]));
              }
          }
      }
  }
  
  D.resize(4*F.rows(), 4*V.rows());
  DA.resize(4*V.rows(), 4*F.rows());
  D.setFromTriplets(IJV_D.begin(), IJV_D.end());
  DA.setFromTriplets(IJV_DA.begin(), IJV_DA.end());

  D.makeCompressed();
  DA.makeCompressed();
}
