// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "triangle_triangle_adjacency.h"
#include "is_edge_manifold.h"
#include "all_edges.h"
#include <map>
#include <algorithm>

template <typename Scalar, typename Index>
IGL_INLINE void igl::triangle_triangle_adjacency_preprocess(const Eigen::PlainObjectBase<Scalar>& /*V*/,
                                   const Eigen::PlainObjectBase<Index>& F,
                                   std::vector<std::vector<int> >& TTT)
{
  for(int f=0;f<F.rows();++f)
    for (int i=0;i<F.cols();++i)
    {
      // v1 v2 f ei
      int v1 = F(f,i);
      int v2 = F(f,(i+1)%F.cols());
      if (v1 > v2) std::swap(v1,v2);
      std::vector<int> r(4);
      r[0] = v1; r[1] = v2;
      r[2] = f;  r[3] = i;
      TTT.push_back(r);
    }
  std::sort(TTT.begin(),TTT.end());
}

// Extract the face adjacencies
template <typename Index>
IGL_INLINE void igl::triangle_triangle_adjacency_extractTT(const Eigen::PlainObjectBase<Index>& F,
                                  std::vector<std::vector<int> >& TTT,
                                  Eigen::PlainObjectBase<Index>& TT)
{
  TT = Eigen::PlainObjectBase<Index>::Constant((int)(F.rows()),F.cols(),-1);

  for(int i=1;i<(int)TTT.size();++i)
  {
    std::vector<int>& r1 = TTT[i-1];
    std::vector<int>& r2 = TTT[i];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
    {
      TT(r1[2],r1[3]) = r2[2];
      TT(r2[2],r2[3]) = r1[2];
    }
  }
}

// Extract the face adjacencies indices (needed for fast traversal)
template <typename Index>
IGL_INLINE void igl::triangle_triangle_adjacency_extractTTi(const Eigen::PlainObjectBase<Index>& F,
                                   std::vector<std::vector<int> >& TTT,
                                   Eigen::PlainObjectBase<Index>& TTi)
{
  TTi = Eigen::PlainObjectBase<Index>::Constant((int)(F.rows()),F.cols(),-1);

  for(int i=1;i<(int)TTT.size();++i)
  {
    std::vector<int>& r1 = TTT[i-1];
    std::vector<int>& r2 = TTT[i];
    if ((r1[0] == r2[0]) && (r1[1] == r2[1]))
    {
      TTi(r1[2],r1[3]) = r2[3];
      TTi(r2[2],r2[3]) = r1[3];
    }
  }
}

// Compute triangle-triangle adjacency
template <typename Scalar, typename Index>
IGL_INLINE void igl::triangle_triangle_adjacency(const Eigen::PlainObjectBase<Scalar>& V,
                        const Eigen::PlainObjectBase<Index>& F,
                        Eigen::PlainObjectBase<Index>& TT)
{
  //assert(igl::is_edge_manifold(V,F));
  std::vector<std::vector<int> > TTT;

  triangle_triangle_adjacency_preprocess(V,F,TTT);
  triangle_triangle_adjacency_extractTT(F,TTT,TT);
}

// Compute triangle-triangle adjacency with indices
template <typename Scalar, typename Index>
IGL_INLINE void igl::triangle_triangle_adjacency(const Eigen::PlainObjectBase<Scalar>& V,
                        const Eigen::PlainObjectBase<Index>& F,
                        Eigen::PlainObjectBase<Index>& TT,
                        Eigen::PlainObjectBase<Index>& TTi)
{
  //assert(igl::is_edge_manifold(V,F));
  std::vector<std::vector<int> > TTT;

  triangle_triangle_adjacency_preprocess(V,F,TTT);
  triangle_triangle_adjacency_extractTT(F,TTT,TT);
  triangle_triangle_adjacency_extractTTi(F,TTT,TTi);
}

template <
  typename DerivedF, 
  typename TTIndex, 
  typename TTiIndex>
  IGL_INLINE void igl::triangle_triangle_adjacency(
    const Eigen::PlainObjectBase<DerivedF> & F,
    std::vector<std::vector<std::vector<TTIndex> > > & TT,
    std::vector<std::vector<std::vector<TTiIndex> > > & TTi)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  assert(F.cols() == 3 && "Faces must be triangles");
  typedef typename DerivedF::Index Index;
  // number of faces
  const int m = F.rows();
  // All occurances of directed edges
  MatrixXi E;
  all_edges(F,E);
  assert(E.rows() == 3*m);
  // uE2E[i] --> {j,k,...} means unique edge i corresponds to face edges j and
  // k (where j-edge comes is the j/m edge of face j%m)
  map<pair<Index,Index>,vector<Index> > uE2E;
  for(int e = 0;e<E.rows();e++)
  {
    Index i = E(e,0);
    Index j = E(e,1);
    if(i<j)
    {
      uE2E[pair<Index,Index>(i,j)].push_back(e);
    }else
    {
      uE2E[pair<Index,Index>(j,i)].push_back(e);
    }
  }
  // E2E[i] --> {j,k,...} means face edge i corresponds to other faces edges j
  // and k
  TT.resize (m,vector<vector<Index> >(F.cols()));
  TTi.resize(m,vector<vector<Index> >(F.cols()));
  for(int e = 0;e<E.rows();e++)
  {
    const Index i = E(e,0);
    const Index j = E(e,1);
    const Index f = e%m;
    const Index c = e/m;
    const vector<Index> & N = 
      i<j ? uE2E[pair<Index,Index>(i,j)] : uE2E[pair<Index,Index>(j,i)];
    for(const auto & ne : N)
    {
      const Index nf = ne%m;
      const Index nc = ne/m;
      TT[f][c].push_back(nf);
      TTi[f][c].push_back(nc);
    }
  }
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::triangle_triangle_adjacency<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
// generated by autoexplicit.sh
template void igl::triangle_triangle_adjacency<Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&);
template void igl::triangle_triangle_adjacency<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#endif
