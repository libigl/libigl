// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "hausdorff.h"
#include "point_mesh_squared_distance.h"
#include "AABB.h"
#include "upsample.h"
#include "slice_mask.h"
#include "remove_unreferenced.h"
#include "EPS.h"

template <
  typename DerivedVA, 
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename Scalar>
IGL_INLINE void igl::hausdorff(
  const Eigen::PlainObjectBase<DerivedVA> & OVA, 
  const Eigen::PlainObjectBase<DerivedFA> & OFA,
  const Eigen::PlainObjectBase<DerivedVB> & VB, 
  const Eigen::PlainObjectBase<DerivedFB> & FB,
  const Scalar eps,
  Scalar & lower,
  Scalar & upper)
{
  using namespace Eigen;
  using namespace igl;
  using namespace std;
  assert(VA.cols() == 3 && "VA should contain 3d points");
  assert(FA.cols() == 3 && "FA should contain triangles");
  assert(VB.cols() == 3 && "VB should contain 3d points");
  assert(FB.cols() == 3 && "FB should contain triangles");
  assert(eps>0 && "eps should be a strictly positive number");
  const auto inf = std::numeric_limits<Scalar>::infinity();
  DerivedVA VA = OVA;
  DerivedFA FA = OFA;
  // Precompute acceleration data structure for B
  AABB<DerivedVB,3> tree;
  tree.init(VB,FB);
  lower = 0;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorXS;
  VectorXS DV;
  Eigen::Matrix<int,Eigen::Dynamic,1> DI;
  Eigen::Matrix<Scalar,Eigen::Dynamic,3> DC;
  while(true)
  {
    // Compute distance to each vertex in VA
    tree.squared_distance(VB,FB,VA,DV,DI,DC);
    DV = DV.array().sqrt();
    // Compute upper bound for each face in FA
    VectorXS U(FA.rows(),1);
    for(int i = 0;i<FA.rows();i++)
    {
      // Compute edge lengths
      Eigen::Matrix<Scalar,3,1> e;
      for(int c = 0;c<3;c++)
      {
        e(c) = (VA.row(FA(i,(c+1)%3)) - VA.row(FA(i,(c+2)%3))).norm();
      }
      // Semiperimeter
      Scalar s = e.array().sum()/2.0;
      // Area
      Scalar A = sqrt(s*(s-e(0))*(s-e(1))*(s-e(2)));
      // Circumradius
      Scalar R = e(0)*e(1)*e(2)/(4.0*A);
      // Inradius
      Scalar r = A/s;
      // Barycenter
      Eigen::Matrix<typename DerivedVA::Scalar,1,3> B = 
        (VA.row(FA(i,0))+ VA.row(FA(i,1))+ VA.row(FA(i,2)))/3.0;
      // These still sometimes contribute
      Scalar u1 = inf;
      Scalar u2 = 0;
      // u3 is a _huge_ improvement over u1, u2
      Scalar u3 = 0;
      // u4 _is contributing_ but it's rather expensive
      Scalar u4 = inf;
      for(int c = 0;c<3;c++)
      {
        const Scalar dic = DV(FA(i,c));
        lower = std::max(lower,dic);
        u1 = std::min(u1,dic+max(e((c+1)%3), e((c+2)%3)));
        u2 = std::max(u2,dic);
        u3 = std::max(u3,dic);
        u3 = std::max(u3, (B-DC.row(FA(i,c))).norm());
        {
          Scalar u4c = 0;
          for(int o = 0;o<3;o++)
          {
            u4c = 
              std::max(u4c,c==o?dic:(VA.row(FA(i,c))-DC.row(FA(i,o))).norm());
          }
          u4 = std::min(u4,u4c);
        }
      }
      u2 += ( s-r > 2.*R ? R : e.maxCoeff()/2.0 );
      U(i) = std::min(std::min(std::min(u1,u2),u3),u4);
    }
    upper = U.maxCoeff();
    //cout<<std::setprecision(17) << FA.rows()<<"\t"<<lower<< "\t"<<upper<<endl;
    if(upper-lower <= eps)
    {
      return;
    }
    // Remove faces with too small upper bounds to impact result
    slice_mask(DerivedFA(FA),(U.array()>lower).eval(),1,FA);
    {
      Eigen::Matrix<typename DerivedFA::Index,Eigen::Dynamic,1> _;
      remove_unreferenced(DerivedVA(VA),DerivedFA(FA),VA,FA,_);
    }
    // Upsample all (remaining) faces in A
    upsample(DerivedVA(VA),DerivedFA(FA),VA,FA);
  }
}

template <
  typename DerivedVA, 
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename Scalar>
IGL_INLINE void igl::hausdorff(
  const Eigen::PlainObjectBase<DerivedVA> & VA, 
  const Eigen::PlainObjectBase<DerivedFA> & FA,
  const Eigen::PlainObjectBase<DerivedVB> & VB, 
  const Eigen::PlainObjectBase<DerivedFB> & FB,
  Scalar & d)
{
  using namespace Eigen;
  auto X = VA.colwise().maxCoeff().array().max(VB.colwise().maxCoeff().array());
  auto N = VA.colwise().minCoeff().array().min(VB.colwise().minCoeff().array());
  const double bbd = (X-N).matrix().norm();
  const double eps = igl::EPS<Scalar>()*bbd;
  double lab,uab,lba,uba;
  igl::hausdorff(VA,FA,VB,FB,eps,lab,uab);
  igl::hausdorff(VB,FB,VA,FA,eps,lba,uba);
  d = std::max(lab,lba);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::hausdorff<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, double>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, double&);
#endif
