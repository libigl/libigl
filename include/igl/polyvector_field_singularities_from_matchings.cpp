// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2016 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include <iostream>
#include <igl/polyvector_field_singularities_from_matchings.h>
#include <igl/is_border_vertex.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/list_to_matrix.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>

template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedTT>
void igl::polyvector_field_one_ring_matchings(const Eigen::PlainObjectBase<DerivedV> &V,
                                              const Eigen::PlainObjectBase<DerivedF> &F,
                                              const std::vector<std::vector<VFType> >& VF,
                                              const Eigen::MatrixXi& E2F,
                                              const Eigen::MatrixXi& F2E,
                                              const Eigen::PlainObjectBase<DerivedTT>& TT,
                                              const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                              const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                              const int vi,
                                              Eigen::MatrixXi &mvi,
                                              Eigen::VectorXi &fi)
{
  int half_degree = match_ab.cols();
  mvi.resize(VF[vi].size()+1,half_degree);
  fi.resize(VF[vi].size()+1,1);
  //start from one face
  //first, check if the vertex is on a boundary
  //then there must be two faces that are on the boundary
  //(other cases not supported)

  int fstart = -1;
  int ind = 0;
  for (int i =0; i<VF[vi].size(); ++i)
  {
    int fi = VF[vi][i];
    for (int  j=0; j<3; ++j)
      if (F(fi,j)==vi && TT(fi,j) == -1)
      {
        ind ++;
        fstart = fi;
        //        break;
      }
  }
  if (ind >1 )
  {
    std::cerr<<"igl::polyvector_field_one_ring_matchings -- vertex "<<vi<< " is on an unusual boundary"<<std::endl;
    exit(1);
  }
  if (fstart == -1)
    fstart = VF[vi][0];
  int current_face = fstart;
  int i =0;
  fi[i] = current_face;
  for (int j=0; j<half_degree; ++j)
    mvi(i,j) = j;

  int next_face = -1;
  while (next_face != fstart && current_face!=-1)
  {
    // look for the vertex
    int j=-1;
    for (unsigned z=0; z<3; ++z)
      if (F(current_face,(z+1)%3) == vi)
      {
        j=z;
        break;
      }
    assert(j!=-1);

    next_face = TT(current_face, j);
    ++i;

    if (next_face == -1)
      mvi.row(i).setConstant(-1);
    else
    {
      // look at the edge between the two faces
      const int &current_edge = F2E(current_face,j);

      for (int k=0; k<half_degree; ++k)
      {
        // check its orientation to determine whether match_ab or match_ba should be used
        if ((E2F(current_edge,0) == current_face) &&
            (E2F(current_edge,1) == next_face) )
        {
          //look at match_ab
          mvi(i,k) = match_ab(current_edge,(mvi(i-1,k))%half_degree);
        }
        else
        {
          assert((E2F(current_edge,1) == current_face) &&
                 (E2F(current_edge,0) == next_face));
          //look at match_ba
          mvi(i,k) = match_ba(current_edge,(mvi(i-1,k))%half_degree);
        }
        if (mvi(i-1,k)>=half_degree)
          mvi(i,k) = (mvi(i,k)+half_degree)%(2*half_degree);
      }
    }
    current_face = next_face;
    fi[i] = current_face;
  }
}

template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedS>
IGL_INLINE void igl::polyvector_field_singularities_from_matchings(
                                                                   const Eigen::PlainObjectBase<DerivedV> &V,
                                                                   const Eigen::PlainObjectBase<DerivedF> &F,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularities)
{

  std::vector<bool> V_border = igl::is_border_vertex(V,F);
  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(F,TT,TTi);

  Eigen::MatrixXi E, E2F, F2E;
  igl::edge_topology(V,F,E,F2E,E2F);

  igl::polyvector_field_singularities_from_matchings(V, F, V_border, VF, TT, E2F, F2E, match_ab, match_ba, singularities);
}

template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedS>
IGL_INLINE void igl::polyvector_field_singularities_from_matchings(
                                                                   const Eigen::PlainObjectBase<DerivedV> &V,
                                                                   const Eigen::PlainObjectBase<DerivedF> &F,
                                                                   const std::vector<bool> &V_border,
                                                                   const std::vector<std::vector<VFType> > &VF,
                                                                   const Eigen::MatrixXi &TT,
                                                                   const Eigen::MatrixXi &E2F,
                                                                   const Eigen::MatrixXi &F2E,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularities)
{

  int numV = V.rows();

  std::vector<int> singularities_v;
  int half_degree = match_ab.cols();
  for (int vi =0; vi<numV; ++vi)
  {
    ///check that is on border..
    if (V_border[vi])
      continue;
    Eigen::VectorXi fi;
    Eigen::MatrixXi mvi;
    igl::polyvector_field_one_ring_matchings(V, F, VF, E2F, F2E, TT, match_ab, match_ba, vi, mvi, fi);

    int num = fi.size();
    //pick one of the vectors to check for singularities
    for (int vector_to_match = 0; vector_to_match < half_degree; ++vector_to_match)
    {
      if(mvi(num-1,vector_to_match) != mvi(0,vector_to_match))
      {
        singularities_v.push_back(vi);
        break;
      }
    }
  }
  std::sort(singularities_v.begin(), singularities_v.end());
  auto last = std::unique(singularities_v.begin(), singularities_v.end());
  singularities_v.erase(last, singularities_v.end());

  igl::list_to_matrix(singularities_v, singularities);
}


template <typename DerivedV, typename DerivedF, typename DerivedM, typename DerivedS>
IGL_INLINE void igl::polyvector_field_singularities_from_matchings(
                                                                   const Eigen::PlainObjectBase<DerivedV> &V,
                                                                   const Eigen::PlainObjectBase<DerivedF> &F,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularity_indices)
{

  std::vector<bool> V_border = igl::is_border_vertex(V,F);
  std::vector<std::vector<int> > VF, VFi;
  igl::vertex_triangle_adjacency(V,F,VF,VFi);

  Eigen::MatrixXi TT, TTi;
  igl::triangle_triangle_adjacency(V,F,TT,TTi);

  Eigen::MatrixXi E, E2F, F2E;
  igl::edge_topology(V,F,E,F2E,E2F);

  igl::polyvector_field_singularities_from_matchings(V, F, V_border, VF, TT, E2F, F2E, match_ab, match_ba, singularities, singularity_indices);
}

template <typename DerivedV, typename DerivedF, typename DerivedM, typename VFType, typename DerivedS>
IGL_INLINE void igl::polyvector_field_singularities_from_matchings(
                                                                   const Eigen::PlainObjectBase<DerivedV> &V,
                                                                   const Eigen::PlainObjectBase<DerivedF> &F,
                                                                   const std::vector<bool> &V_border,
                                                                   const std::vector<std::vector<VFType> > &VF,
                                                                   const Eigen::MatrixXi &TT,
                                                                   const Eigen::MatrixXi &E2F,
                                                                   const Eigen::MatrixXi &F2E,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ab,
                                                                   const Eigen::PlainObjectBase<DerivedM> &match_ba,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularities,
                                                                   Eigen::PlainObjectBase<DerivedS> &singularity_indices)
{
  igl::polyvector_field_singularities_from_matchings(V, F, V_border, VF, TT, E2F, F2E, match_ab, match_ba, singularities);

  singularity_indices.setZero(singularities.size(), 1);

  //get index from first vector only
  int vector_to_match = 0;
  for (int i =0; i<singularities.size(); ++i)
  {
    int vi = singularities[i];

    // Eigen::VectorXi mvi,fi;
    // igl::polyvector_field_one_ring_matchings(V, F, VF, E2F, F2E, TT, match_ab, match_ba, vi, vector_to_match, mvi, fi);
    Eigen::VectorXi fi;
    Eigen::MatrixXi mvi;
    igl::polyvector_field_one_ring_matchings(V, F, VF, E2F, F2E, TT, match_ab, match_ba, vi, mvi, fi);

    singularity_indices[i] = (mvi(mvi.rows()-1,vector_to_match) - vector_to_match);
  }

}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::polyvector_field_singularities_from_matchings<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
