// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <igl/polyvector_field_matchings.h>
#include <igl/edge_topology.h>
#include <igl/per_face_normals.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/rotation_matrix_from_directions.h>

template <typename DerivedS, typename DerivedV, typename DerivedM>
IGL_INLINE void igl::polyvector_field_matching(
                                               const Eigen::PlainObjectBase<DerivedS>& _ua,
                                               const Eigen::PlainObjectBase<DerivedS>& _ub,
                                               const Eigen::PlainObjectBase<DerivedV>& na,
                                               const Eigen::PlainObjectBase<DerivedV>& nb,
                                               const Eigen::PlainObjectBase<DerivedV>& e,
                                               bool match_with_curl,
                                               Eigen::PlainObjectBase<DerivedM>& mab,
                                               Eigen::PlainObjectBase<DerivedM>& mba,
                                               bool is_symmetric)
{
  // make sure the matching preserve ccw order of the vectors across the edge
  // 1) order vectors in a, ccw  e.g. (0,1,2,3)_a not ccw --> (0,3,2,1)_a ccw
  // 2) order vectors in b, ccw  e.g. (0,1,2,3)_b not ccw --> (0,2,1,3)_b ccw
  // 3) the vectors in b that match the ordered vectors in a (in this case  (0,3,2,1)_a ) must be a circular shift of the ccw ordered vectors in b  - so we need to explicitly check only these matchings to find the best ccw one, there are N of them
  int hN = _ua.cols()/3;
  int N;
  if (is_symmetric)
    N = 2*hN;
  else
    N = hN;

  Eigen::Matrix<typename DerivedS::Scalar,1,Eigen::Dynamic> ua (1,N*3);
  Eigen::Matrix<typename DerivedS::Scalar,1,Eigen::Dynamic> ub (1,N*3);
  if (is_symmetric)
  {
    ua <<_ua, -_ua;
    ub <<_ub, -_ub;
  }
  else
  {
    ua =_ua;
    ub =_ub;
  }

  Eigen::Matrix<typename DerivedM::Scalar,Eigen::Dynamic,1> order_a, order_b;
  igl::sort_vectors_ccw(ua, na, order_a);
  igl::sort_vectors_ccw(ub, nb, order_b);

  //checking all possible circshifts of order_b as matches for order_a
  Eigen::Matrix<typename DerivedM::Scalar,Eigen::Dynamic,Eigen::Dynamic> all_matches(N,N);
  Eigen::Matrix<typename DerivedS::Scalar,1,Eigen::Dynamic> all_scores(1,N);
  for (int s =0; s< N; ++s)
  {
    all_matches.row(s) = order_b;
    typename DerivedS::Scalar current_score=0;
    for (int i=0; i< N; ++i)
    {
      if (match_with_curl)
        current_score += fabs(ua.segment(order_a[i]*3, 3).dot(e) - ub.segment(order_b[i]*3, 3).dot(e));
      else
      {
        Eigen::Matrix<typename DerivedS::Scalar,3,1> na_ = na.transpose();
        Eigen::Matrix<typename DerivedS::Scalar,3,1> nb_ = nb.transpose();

        Eigen::Matrix<typename DerivedS::Scalar,1,3> uaRot = igl::rotation_matrix_from_directions(na_, nb_)*ua.segment(order_a[i]*3, 3).transpose();
        current_score += (uaRot-ub.segment(order_b[i]*3, 3)).norm();
      }
    }
    all_scores[s] = current_score;
    // do a circshift on order_b to check the next preserving matching
    int temp = order_b[0];
    for (int i =0; i< N-1; ++i)
      order_b[i] = order_b[i+1];
    order_b(N-1) = temp;
  }
  Eigen::Matrix<typename DerivedM::Scalar,1,Eigen::Dynamic> best_matching_for_sorted_a;
  int best_i;
  all_scores.minCoeff(&best_i);
  best_matching_for_sorted_a = all_matches.row(best_i);
  // best_matching_for_sorted_a is the matching for the sorted vectors in a
  // get the matching for the initial (unsorted) vectors
  mab.resize(1,N);
  for (int i=0; i< N; ++i)
    mab[order_a[i]] = best_matching_for_sorted_a[i];

  //mab contains the best matching a->b, get the opposite too
  mba.resize(1, N);
  for (int i=0; i< N; ++i)
    mba[mab[i]] = i;

  if (is_symmetric)
  {
  mab = mab.head(hN);
  mba = mba.head(hN);
}
}


template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedM>
IGL_INLINE void igl::polyvector_field_matchings(
                                                                     const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                                     const Eigen::PlainObjectBase<DerivedV>&V,
                                                                     const Eigen::PlainObjectBase<DerivedF>&F,
                                                                     const Eigen::PlainObjectBase<DerivedE>&E,
                                                                     const Eigen::PlainObjectBase<DerivedV>& FN,
                                                                     const Eigen::MatrixXi &E2F,
                                                                     bool match_with_curl,
                                                bool is_symmetric,
                                                                     Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                Eigen::PlainObjectBase<DerivedM>& match_ba)
{
  int numEdges = E.rows();
  int half_degree = sol3D.cols()/3;

  Eigen::VectorXi isBorderEdge;
  isBorderEdge.setZero(numEdges,1);
  for(unsigned i=0; i<numEdges; ++i)
  {
    if ((E2F(i,0) == -1) || ((E2F(i,1) == -1)))
      isBorderEdge[i] = 1;
  }

  match_ab.setZero(numEdges, half_degree);
  match_ba.setZero(numEdges, half_degree);

  for (int ei=0; ei<numEdges; ++ei)
  {
    if (isBorderEdge[ei])
      continue;
    // the two faces of the flap
    int a = E2F(ei,0);
    int b = E2F(ei,1);

    Eigen::Matrix<typename DerivedV::Scalar, 1, Eigen::Dynamic> ce = (V.row(E(ei,1))-V.row(E(ei,0))).normalized().template cast<typename DerivedV::Scalar>();

    Eigen::Matrix<typename DerivedM::Scalar, 1, Eigen::Dynamic> mab, mba;
    igl::polyvector_field_matching(sol3D.row(a).eval(),
                                   sol3D.row(b).eval(),
                                   FN.row(a).eval(),
                                   FN.row(b).eval(),
                                   ce,
                                   match_with_curl,
                                   mab,
                                   mba,
                                   is_symmetric);

    match_ab.row(ei) = mab;
    match_ba.row(ei) = mba;
  }
}


template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedE, typename DerivedM, typename DerivedC>
IGL_INLINE typename DerivedC::Scalar igl::polyvector_field_matchings(
                                                                     const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                                     const Eigen::PlainObjectBase<DerivedV>&V,
                                                                     const Eigen::PlainObjectBase<DerivedF>&F,
                                                                     const Eigen::PlainObjectBase<DerivedE>&E,
                                                                     const Eigen::PlainObjectBase<DerivedV>& FN,
                                                                     const Eigen::MatrixXi &E2F,
                                                                     bool match_with_curl,
                                                                     bool is_symmetric,
                                                                     Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                                     Eigen::PlainObjectBase<DerivedM>& match_ba,
                                                                     Eigen::PlainObjectBase<DerivedC>& curl)
{
  int numEdges = E.rows();
  int half_degree = sol3D.cols()/3;

  Eigen::VectorXi isBorderEdge;
  isBorderEdge.setZero(numEdges,1);
  for(unsigned i=0; i<numEdges; ++i)
  {
    if ((E2F(i,0) == -1) || ((E2F(i,1) == -1)))
      isBorderEdge[i] = 1;
  }

  igl::polyvector_field_matchings(sol3D, V, F, E, FN, E2F, match_with_curl, is_symmetric, match_ab, match_ba);
  curl.setZero(numEdges,1);
  typename DerivedC::Scalar meanCurl = 0;
  for (int ei=0; ei<numEdges; ++ei)
  {
    if (isBorderEdge[ei])
      continue;
    // the two faces of the flap
    int a = E2F(ei,0);
    int b = E2F(ei,1);

    const Eigen::Matrix<typename DerivedM::Scalar, 1, Eigen::Dynamic> &mab = match_ab.row(ei);
    const Eigen::Matrix<typename DerivedM::Scalar, 1, Eigen::Dynamic> &mba = match_ba.row(ei);

    Eigen::Matrix<typename DerivedS::Scalar, 1, Eigen::Dynamic> matched;
    matched.resize(1, 3*half_degree);
    for (int i = 0; i<half_degree; ++i)
    {
      int sign = 1;
      int m = mab[i];
      if (is_symmetric)
      {
        sign = (mab[i]<half_degree)?1:-1;
        m = m%half_degree;
    }
      matched.segment(i*3, 3) = sign*sol3D.row(b).segment(m*3, 3);
    }

    Eigen::Matrix<typename DerivedV::Scalar, 1, Eigen::Dynamic> ce = (V.row(E(ei,1))-V.row(E(ei,0))).normalized().template cast<typename DerivedV::Scalar>();
    typename DerivedC::Scalar avgCurl = 0;
    for (int i = 0; i<half_degree; ++i)
      avgCurl += fabs(sol3D.row(a).segment(i*3, 3).dot(ce) - matched.segment(i*3, 3).dot(ce));

    avgCurl = avgCurl/half_degree;


    curl[ ei ] = avgCurl;

    meanCurl+= avgCurl;

  }

  meanCurl /= 1.*(numEdges - isBorderEdge.sum());
  return meanCurl;
}

template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedC>
IGL_INLINE typename DerivedC::Scalar igl::polyvector_field_matchings(
                                                                     const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                                     const Eigen::PlainObjectBase<DerivedV>&V,
                                                                     const Eigen::PlainObjectBase<DerivedF>&F,
                                                                     bool match_with_curl,
                                                                     bool is_symmetric,
                                                                     Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                                     Eigen::PlainObjectBase<DerivedM>& match_ba,
                                                                     Eigen::PlainObjectBase<DerivedC>& curl)
{
  Eigen::MatrixXi E, E2F, F2E;
  igl::edge_topology(V,F,E,F2E,E2F);

//#warning "Poor templating of igl::polyvector_field_matchings forces FN to be DerivedV (this could cause issues if DerivedV has fixed number of rows)"
  DerivedV FN;
  igl::per_face_normals(V,F,FN);

  return igl::polyvector_field_matchings(sol3D, V, F, E, FN, E2F, match_with_curl, is_symmetric, match_ab, match_ba, curl);
}

template <typename DerivedS, typename DerivedV, typename DerivedF, typename DerivedM>
IGL_INLINE void igl::polyvector_field_matchings(
                                                const Eigen::PlainObjectBase<DerivedS>& sol3D,
                                                const Eigen::PlainObjectBase<DerivedV>&V,
                                                const Eigen::PlainObjectBase<DerivedF>&F,
                                                bool match_with_curl,
                                                bool is_symmetric,
                                                Eigen::PlainObjectBase<DerivedM>& match_ab,
                                                Eigen::PlainObjectBase<DerivedM>& match_ba)
{
  Eigen::MatrixXi E, E2F, F2E;
  igl::edge_topology(V,F,E,F2E,E2F);

  DerivedV FN;
  igl::per_face_normals(V,F,FN);

  igl::polyvector_field_matchings(sol3D, V, F, E, FN, E2F, match_with_curl, is_symmetric, match_ab, match_ba);
}


#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::polyvector_field_matchings<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template Eigen::Matrix<double, -1, 1, 0, -1, 1>::Scalar igl::polyvector_field_matchings<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, bool, bool, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1> >&);
#endif
