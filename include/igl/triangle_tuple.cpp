// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <zhongshi@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "triangle_tuple.h"

namespace igl
{
  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE void triangle_tuple_switch_vert(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    along = !along;
  };

  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE void triangle_tuple_switch_edge(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    ei = (ei + ((!along) ? 1 : 2)) % 3;
    triangle_tuple_switch_vert(fi, ei, along, F, FF, FFi);
  };

  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE void triangle_tuple_switch_face(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    int fin = FF(fi, ei);
    int ein = FFi(fi, ei);

    fi = fin;
    ei = ein;
    triangle_tuple_switch_vert(fi, ei, along, F, FF, FFi);
  };

  template<typename DerivedF, typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_vert(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    assert(fi <= F.rows());
    assert(fi >= 0);
    assert(ei >= 0);
    assert(ei <= 2);

    // legacy edge indexing
    return F(fi, (!along) ? (ei + 1) % 3 : ei);
  };

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_edge(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    return ei;
  };

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE int triangle_tuple_get_face(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    return fi;
  };

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE bool triangle_tuple_next_in_one_ring(
      int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    if(triangle_tuple_is_on_boundary(fi, ei, along, F, FF, FFi))
    {
      do
      {
        triangle_tuple_switch_face(fi, ei, along, F, FF, FFi);
        triangle_tuple_switch_edge(fi, ei, along, F, FF, FFi);
      } while(!triangle_tuple_is_on_boundary(fi, ei, along, F, FF, FFi));
      triangle_tuple_switch_edge(fi, ei, along, F, FF, FFi);
      return false;
    }
    else
    {
      triangle_tuple_switch_face(fi, ei, along, F, FF, FFi);
      triangle_tuple_switch_edge(fi, ei, along, F, FF, FFi);
      return true;
    }
  };

  template<typename DerivedF,
           typename DerivedFF,
           typename DerivedFFi>
  IGL_INLINE bool triangle_tuple_is_on_boundary(
      const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedF> &F,
      const Eigen::MatrixBase<DerivedFF> &FF,
      const Eigen::MatrixBase<DerivedFFi> &FFi)
  {
    return F(fi, ei) == -1;
  };

  IGL_INLINE bool triangle_tuples_equal(
      const int &f1, const int &e1, const bool &a1,
      const int &f2, const int &e2, const bool &a2)
  {
    return f1 == f2 &&
           e1 == e2 &&
           a1 == a2;
  };
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::triangle_tuple_get_vert<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int const&, int const&, bool const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template bool igl::triangle_tuple_next_in_one_ring<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int &, int &, bool &, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
#endif
