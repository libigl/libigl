// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Zhongshi Jiang <zhongshi@cims.nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "tetrahedron_tuple.h"
#include "tetrahedron_tetrahedron_adjacency.h"
namespace igl
{
  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tet_tuple_switch_vert(
      int &ti, int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    along = !along;
  };

  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tet_tuple_switch_edge(
      int &ti, int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    ei = (ei + (along ? 2 : 1)) % 3;
    tet_tuple_switch_vert(ti, fi, ei, along, T, TT, TTif, TTie);
  };

  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tet_tuple_switch_face(
      int &ti, int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    fi = tetrahedron_local_FF(fi, (ei + 2) % 3);
    tet_tuple_switch_vert(ti, fi, ei, along, T, TT, TTif, TTie);
  };

  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE void tet_tuple_switch_tet(
      int &ti, int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    if(tet_tuple_is_on_boundary(ti, fi, ei, along, T, TT, TTif, TTie)) return;
    int tin = TT(ti, fi);
    int fin = TTif(ti, fi);
    int ein = TTie(ti, 3 * fi + ei);

    ti = tin;
    fi = fin;
    ei = ein;
    tet_tuple_switch_vert(ti, fi, ei, along, T, TT, TTif, TTie);
  };

  template<typename DerivedT, typename DerivedTT,
           typename DerivedTTif, typename DerivedTTie>
  IGL_INLINE int tet_tuple_get_vert(
      const int &ti, const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    assert(ti >= 0);
    assert(ti < T.rows());
    assert(fi <= 3);
    assert(fi >= 0);
    assert(ei >= 0);
    assert(ei <= 2);

    // legacy edge indexing
    return T(ti,
             (tetrahedron_local_FF)(fi,
                                    along ? ei: (ei + 1) % 3));
  };

  template<typename DerivedT,
           typename DerivedTT,
           typename DerivedTTif,
           typename DerivedTTie>
  IGL_INLINE int tet_tuple_get_edge(
      const int &ti, const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    return ei;
  };

  template<typename DerivedT,
           typename DerivedTT,
           typename DerivedTTif,
           typename DerivedTTie>
  IGL_INLINE int tet_tuple_get_face(
      const int &ti, const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    return fi;
  };

  template<typename DerivedT,
           typename DerivedTT,
           typename DerivedTTif,
           typename DerivedTTie>
  IGL_INLINE int tet_tuple_get_tet(
      const int &ti, const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    return ti;
  };

  template<typename DerivedT,
           typename DerivedTT,
           typename DerivedTTif,
           typename DerivedTTie>
  IGL_INLINE bool tet_tuple_next_in_one_ring(
      int &ti, int &fi, int &ei, bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    if(tet_tuple_is_on_boundary(ti, fi, ei, along, T, TT, TTif, TTie))
    {
      do
      {
        tet_tuple_switch_face(ti, fi, ei, along, T, TT, TTif, TTie);
        tet_tuple_switch_tet(ti, fi, ei, along, T, TT, TTif, TTie);
        tet_tuple_switch_face(ti, fi, ei, along, T, TT, TTif, TTie);
        tet_tuple_switch_edge(ti, fi, ei, along, T, TT, TTif, TTie);
      } while(!tet_tuple_is_on_boundary(ti, fi, ei, along, T, TT, TTif,
                                        TTie));
      tet_tuple_switch_edge(ti, fi, ei, along, T, TT, TTif, TTie);
      return false;
    }
    else
    {
      tet_tuple_switch_face(ti, fi, ei, along, T, TT, TTif, TTie);
      tet_tuple_switch_tet(ti, fi, ei, along, T, TT, TTif, TTie);
      tet_tuple_switch_face(ti, fi, ei, along, T, TT, TTif, TTie);
      tet_tuple_switch_edge(ti, fi, ei, along, T, TT, TTif, TTie);
      return true;
    }
  };

  template<typename DerivedT,
           typename DerivedTT,
           typename DerivedTTif,
           typename DerivedTTie>
  IGL_INLINE bool tet_tuple_is_on_boundary(
      const int &ti, const int &fi, const int &ei, const bool &along,
      const Eigen::MatrixBase<DerivedT> &T,
      const Eigen::MatrixBase<DerivedTT> &TT,
      const Eigen::MatrixBase<DerivedTTif> &TTif,
      const Eigen::MatrixBase<DerivedTTie> &TTie)
  {
    return TT(ti, fi) == -1;
  };

  IGL_INLINE bool tet_tuples_equal(
      const int &t1, const int &f1, const int &e1, const bool &rev1,
      const int &t2, const int &f2, const int &e2, const bool &rev2)
  {
    return t1 == t2 &&
           f1 == f2 &&
           e1 == e2 &&
           rev1 == rev2;
  };
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template int igl::tet_tuple_get_tet<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(const int&, const int&, const int&, const bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template int igl::tet_tuple_get_edge<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(const int&, const int&, const int&, const bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template int igl::tet_tuple_get_face<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(const int&, const int&, const int&, const bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template int igl::tet_tuple_get_vert<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(const int&, const int&, const int&, const bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::tet_tuple_switch_tet<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int&, int&, int&, bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::tet_tuple_switch_edge<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int&, int&, int&, bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::tet_tuple_switch_face<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int&, int&, int&, bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
template void igl::tet_tuple_switch_vert<Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int&, int&, int&, bool&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&);
#endif
