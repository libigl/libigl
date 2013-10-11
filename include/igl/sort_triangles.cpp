#include "sort_triangles.h"
#include "barycenter.h"
#include "sort.h"
#include "sortrows.h"
#include "slice.h"
#include "round.h"
#include "colon.h"
#include "matlab_format.h"
#include "OpenGL_convenience.h"

#include <iostream>

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedMV,
  typename DerivedP,
  typename DerivedFF,
  typename DerivedI>
IGL_INLINE void igl::sort_triangles(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedMV> & MV,
  const Eigen::PlainObjectBase<DerivedP> & P,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;


  // Barycenter, centroid
  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,1> D,sD;
  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,4> BC,PBC;
  barycenter(V,F,BC);
  D = BC*(MV.transpose()*P.transpose().eval().col(2));
  sort(D,1,false,sD,I);

  //// Closest corner
  //Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,1> D,sD;
  //D.setConstant(F.rows(),1,-1e26);
  //for(int c = 0;c<3;c++)
  //{
  //  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,4> C;
  //  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,1> DC;
  //  C.resize(F.rows(),4);
  //  for(int f = 0;f<F.rows();f++)
  //  {
  //    C(f,0) = V(F(f,c),0);
  //    C(f,1) = V(F(f,c),1);
  //    C(f,2) = V(F(f,c),2);
  //    C(f,3) = 1;
  //  }
  //  DC = C*(MV.transpose()*P.transpose().eval().col(2));
  //  D = (DC.array()>D.array()).select(DC,D).eval();
  //}
  //sort(D,1,false,sD,I);

  //// Closest corner with tie breaks
  //Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,3> D,sD,ssD;
  //D.resize(F.rows(),3);
  //for(int c = 0;c<3;c++)
  //{
  //  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,4> C;
  //  C.resize(F.rows(),4);
  //  for(int f = 0;f<F.rows();f++)
  //  {
  //    C(f,0) = V(F(f,c),0);
  //    C(f,1) = V(F(f,c),1);
  //    C(f,2) = V(F(f,c),2);
  //    C(f,3) = 1;
  //  }
  //  D.col(c) = C*(MV.transpose()*P.transpose().eval().col(2));
  //}
  //VectorXi _;
  //sort(D,2,false,sD,_);
  //sortrows(sD,false,ssD,I);


  slice(F,I,1,FF);
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedFF,
  typename DerivedI>
void igl::sort_triangles(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  // Put model, projection, and viewport matrices into double arrays
  Matrix4d MV;
  Matrix4d P;
  glGetDoublev(GL_MODELVIEW_MATRIX,  MV.data());
  glGetDoublev(GL_PROJECTION_MATRIX, P.data());
  if(V.cols() == 3)
  {
    Matrix<typename DerivedV::Scalar, DerivedV::RowsAtCompileTime,4> hV;
    hV.resize(V.rows(),4);
    hV.block(0,0,V.rows(),V.cols()) = V;
    hV.col(3).setConstant(1);
    return sort_triangles(hV,F,MV,P,FF,I);
  }else
  {
    return sort_triangles(V,F,MV,P,FF,I);
  }
}


#include "project.h"
#include <iostream>
template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedFF,
  typename DerivedI>
void igl::sort_triangles_slow(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  Eigen::PlainObjectBase<DerivedFF> & FF,
  Eigen::PlainObjectBase<DerivedI> & I)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  // Barycenter, centroid
  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,1> D,sD;
  Eigen::Matrix<typename DerivedV::Scalar, DerivedF::RowsAtCompileTime,3> BC;
  D.resize(F.rows(),3);
  barycenter(V,F,BC);
  for(int f = 0;f<F.rows();f++)
  {
    Eigen::Matrix<typename DerivedV::Scalar, 3,1> bc,pbc;
    bc = BC.row(f);
    project(bc,pbc);
    D(f) = pbc(2);
  }
  sort(D,1,false,sD,I);
  slice(F,I,1,FF);
}

#ifndef IGL_HEADER_ONLY
// Explicit template instanciation
template void igl::sort_triangles<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::sort_triangles_slow<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
