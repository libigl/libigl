// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "outer_facet.h"
#include "sort.h"
#include "vertex_triangle_adjacency.h"
#include <iostream>
#define IGL_OUTER_FACET_DEBUG

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedN,
  typename DerivedI,
  typename f_type>
IGL_INLINE void igl::outer_facet(
  const Eigen::PlainObjectBase<DerivedV> & V,
  const Eigen::PlainObjectBase<DerivedF> & F,
  const Eigen::PlainObjectBase<DerivedN> & N,
  const Eigen::PlainObjectBase<DerivedI> & I,
  f_type & max_f,
  bool & flip)
{
  using namespace std;
  typedef typename DerivedV::Scalar Scalar;
  typedef typename DerivedV::Index Index;
  typedef 
    typename Eigen::Matrix<Index, DerivedV::RowsAtCompileTime,1> VectorXI;
  typedef 
    typename Eigen::Matrix<Scalar, DerivedV::RowsAtCompileTime,1> VectorXS;
  // "Direct repair of self-intersecting meshes" [Attene 14]
  const Index mi = I.size();
  assert(V.cols() == 3);
  assert(N.cols() == 3);
  Index max_v = -1;
  auto generic_fabs = [&](Scalar val) { 
      if (val >= 0) return val;
      else return -val;
  };
  for(size_t d = 0;d<(size_t)V.cols();d++)
  {
    Scalar max_d = -1e26;
    Scalar max_nd = -1e26;
    for(Index i = 0;i<mi;i++)
    {
      const Index f = I(i);
      const Scalar nd = N(f,d);
      if(generic_fabs(nd)>0)
      {
        for(Index c = 0;c<3;c++)
        {
          const Index v = F(f,c);
          if(v == max_v)
          {
            if(generic_fabs(nd) > max_nd)
            {
              // Just update max face and normal
              max_f = f;
              max_nd = generic_fabs(nd);
              flip = nd<0;
            } else if (generic_fabs(nd) == max_nd) {
                if (nd == max_nd) {
                    if (flip) {
                        max_f = f;
                        max_nd = nd;
                        flip = false;
                    } else if (f > (Index)max_f){
                        max_f = f;
                        max_nd = nd;
                        flip = false;
                    }
                } else {
                    if (flip && f < (Index)max_f) {
                        max_f = f;
                        max_nd = generic_fabs(nd);
                        flip = true;
                    }
                }
            }
          }else
          {
            const Scalar vd = V(v,d);
            if(vd>max_d)
            {
              // update max vertex, face and normal
              max_v = v;
              max_d = vd;
              max_f = f;
              max_nd = generic_fabs(nd);
              flip = nd<0;
            }
          }
        }
      }
    }
    if(max_v >= 0)
    {
      break;
    }
    // if we get here and max_v is still -1 then there were no faces with
    // |N(d)| > 0
  }
#ifdef IGL_OUTER_FACET_DEBUG
  if(max_v <0)
  {
    cerr<<"Very degenerate case, no suitable face found."<<endl;
  }
#endif
  assert(max_v >=0 && "Very degenerate case, no suitable face found.");
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::outer_facet<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::outer_facet<Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 1, -1, -1>, Eigen::Matrix<double, -1, -1, 1, -1, -1>, Eigen::Matrix<int, -1, -1, 1, -1, -1>, unsigned long>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 1, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 1, -1, -1> > const&, unsigned long&, bool&);
template void igl::outer_facet<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
template void igl::outer_facet<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, int>(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int&, bool&);
#endif
