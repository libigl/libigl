// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2022 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "split_nonmanifold.h"
#include "unique_edge_map.h"
#include "connected_components.h"
#include "triangle_triangle_adjacency.h"
#include <cassert>
#include <type_traits>

#include "is_vertex_manifold.h"
#include "matlab_format.h"
#include <iostream>

template <
  typename DerivedF,
  typename DerivedSF,
  typename DerivedSVI
  >
IGL_INLINE void igl::split_nonmanifold(
  const Eigen::MatrixBase<DerivedF> & F,
  Eigen::PlainObjectBase <DerivedSF> & SF,
  Eigen::PlainObjectBase <DerivedSVI> & SVI)
{
  const bool enforce_orientability = true;
#warning "Another parameter whether to try to weld together orientable cut-boundarieS"
  using Scalar = typename DerivedSF::Scalar;
  // Scalar must allow negative values
  static_assert(std::is_signed<Scalar>::value,"Scalar must be signed");
  using MatrixX2I = Eigen::Matrix<Scalar,Eigen::Dynamic,2>;
  using MatrixX3I = Eigen::Matrix<Scalar,Eigen::Dynamic,3>;
  using VectorXI = Eigen::Matrix< Scalar,Eigen::Dynamic,1>;
  MatrixX2I E,uE;
  VectorXI EMAP,uEC,uEE;
  igl::unique_edge_map(F,E,uE,EMAP,uEC,uEE);

  // Mesh as if all edges got cut. 
  MatrixX3I CF = VectorXI::LinSpaced(F.size(),0,F.size()-1).reshaped(F.rows(),F.cols());

  // Every edge starts as a boundary edge
  const int m = F.rows();
  // augh, this is following gptoolbox ordering not triangle_triangle_adjacency
  // ordering.
  MatrixX3I CTF = MatrixX3I::Constant(m,3,-1);
  MatrixX3I CTI = MatrixX3I::Constant(m,3,-1);

  // Is maintaining triangle adjacency enough? Or do we need edge-flaps?
  // After we zip an edge-pair we need to know if there are other unzipped edges
  // adjacent to the zipped edge.
  // Triangle adjacency seems like enough to circulate around the zipped edge
  // vertices and find incident boundary edges, which are candidates for
  // zipping. From a candidate boundary edge we can get to EMAP, which should be
  // enough to find candidate zip-pairs using uE*.

  //// Empty edge-flap data
  //MatrixX2I CEF= MatrixX2I::Constant(E.rows(),2,-1);
  //MatrixX2I CEI= MatrixX2I::Constant(E.rows(),2,-1);

  //std::cout<<igl::matlab_format_index(CF,"CF")<<std::endl;

  // By cutting _all_ edges we also handle non-manifold vertices. We could avoid
  // cutting edges that are both manifold-edges and not incident on a
  // non-manifold vertex. It's hard to imagine really avoiding O(E) work in
  // total though.

  std::vector<bool> seen(E.rows());

  std::vector<std::pair<int,int> > R;

  const auto zip = 
    [&E,&m,&seen,&F,&CF,&R,&CTF,&CTI](const int e1, const int e2, const bool verbose = true)
  {
    assert(!seen[e1]);
    assert(!seen[e2]);
    assert(E(e1,0) == E(e2,1));
    assert(E(e1,1) == E(e2,0));
    const int vs1 = E(e1,0);
    const int vs2 = E(e2,1);
    const int vd1 = E(e1,1);
    const int vd2 = E(e2,0);

    const int f1 = e1%m;
    const int f2 = e2%m;
    const int i1 = e1/m;
    const int i2 = e2/m;

    const int cs1 = CF(f1,(i1+1)%3);
    const int cs2 = CF(f2,(i2+2)%3);
    const int cd1 = CF(f1,(i1+2)%3);
    const int cd2 = CF(f2,(i2+1)%3);
    if(verbose)
    {
    //printf("%d,%d → %d,%d\n",vs1+1,vd1+1,vs2+1,vd2+1);
    //printf("  %d,%d → %d,%d\n",cs1+1,cd1+1,cs2+1,cd2+1);
    printf("  %d → %d\n",f1+1,f2+1);
    }

    const int cs = std::min(cs1,cs2);
    const int cd = std::min(cd1,cd2);
    //printf("  %d,%d → %d,%d\n",cs+1,cd+1,cs+1,cd+1);

    // Record equivalences
    R.emplace_back(cs1,cs);
    R.emplace_back(cd1,cd);
    R.emplace_back(cs2,cs);
    R.emplace_back(cd2,cd);

    CTF(f1,i1) = f2;
    CTI(f1,i1) = i2;
    CTF(f2,i2) = f1;
    CTI(f2,i2) = i1;

    seen[e1] = true;
    seen[e2] = true;

    //std::cout<<igl::matlab_format_index(CF,"CF")<<std::endl;
    //printf("\n");
    
  };


  // consider all unique edges
  for(int u = 0;u<uE.rows();u++)
  {
    const int num_incident = uEC(u+1)-uEC(u);
    assert(num_incident > 0);
    // First edge
    int e1 = uEE(uEC(u));
    // mark boundary edges as seen
    if(num_incident == 1)
    {
      seen[e1] = true;
      continue;
    }
    if(num_incident == 2)
    {
      int e2 = uEE(uEC(u)+1);
      // Check that e2 has opposite orientation of e1
      if(enforce_orientability && E(e1,0) == E(e2,1))
      {
        assert(E(e1,1) == E(e2,0));
        zip(e1,e2,false);
#warning "is `seen` every used?"
        continue;
      }
    }// else. non-manifold. Skip for now
  }

  // Cutting and Stitching: Converting Sets of Polygons to Manifold Surfaces

  // Outer loop over all unique edges.
  for(int u = 0;u<uE.rows();u++)
  {
    const int num_incident = uEC(u+1)-uEC(u);
    assert(num_incident > 0);
    // boundary or internal edge handled already.
    if(num_incident <= 2) { continue; }
    for(int j = uEC(u);j<uEC(u+1);j++)
    {
      const int e1 = uEE(j); // i = j-uEC(u);
      if(seen[e1]) { continue; }
      // We haven't seen this half-edge yet.
      // Try to zip it with an opposite half-edge (and start off a cascade of
      // zipping).
      for(int k = j+1;k<uEC(u+1);k++)
      {
        const int e2 = uEE(k);
        // Can't zip if alread seen
        if(seen[e2]) { continue; }
        if(enforce_orientability && E(e1,0) == E(e2,1))
        {
          assert(E(e1,1) == E(e2,0));
          zip(e1,e2);
          // Don't try any more zips with e1
          break;
        }
      }
    }
  }
  

  // Now if we resolve all the equivalence records we'd have a mesh CF which is
  // "cut" along all original non-manifold edges (and non-manifold vertices
  // separated). This is likely still too aggressive. Every non-manifold edge
  // with 3 incident faces is cut into 3 non-neighboring faces. There's always
  // two of those (with consistent orientation) that can be welded together.

  // We could continue to collect records which conduct these welds along
  // non-manifold edges. We should prioritze welds at "crack" tips: edges whose
  // incident corners have records.
  //

  // We should prioritize welds along chains of non-manifold edges

//#error "Huh? Should we just have cut at non-manifold edges and then considered the incident of connecteced-compontents and non-manifold-edge chains?"
//  // Processing the chains one by one, if a chain has 
//  // if a chain has incident components A,B,C
//  // and A,C are consistently oriented to the chain then we can weld them. 
//
//  A single component could have a shared non-manifold edge with 3 or 4
//  "strips" of faces (like two twisted loop strips out of a shared central
//  mesh colliding through the same non-manifold edge). This has a "manifold"
//  explaination as a single component, but won't be easily merged


  // I couldn't think of a way to do this on the fly. We've collected all the
  // corners that should be mapped to the same vertex as equivalence records
  // (i,j) and we resolve them now using component analysis.
  {
    Eigen::SparseMatrix<bool> A(F.size(),F.size());
    std::vector<Eigen::Triplet<bool> > IJV;
    // Diagonal
    for(int i = 0;i<SF.size();i++) { IJV.emplace_back(i,i,true); }
    // Off-diagonal for each record
    for(const auto & r : R)
    {
      IJV.emplace_back(r.first,r.second,true);
      IJV.emplace_back(r.second,r.first,true);
    }
    // connected compontents
    A.setFromTriplets(IJV.begin(),IJV.end());
    Eigen::VectorXi C,K;
    igl::connected_components(A,C,K);
    for(int i = 0;i<CF.size();i++)
    {
      CF(i) = C(CF(i));
    }
  }

  //{
  //  // Did we correctly maintain CTF and CTI?
  //  MatrixX3I CTF2,CTI2;
  //  std::cout<<igl::matlab_format_index(CF,"CF")<<std::endl;
  //  igl::triangle_triangle_adjacency(CF,CTF2,CTI2);

  //  {
  //    Eigen::PermutationMatrix<3,3> perm(3);
  //    perm.indices() = Eigen::Vector3i(1,2,0);
  //    CTF2 = (CTF2*perm).eval();
  //    CTI2 = (CTI2*perm).eval();
  //    for(int i=0;i<CTI2.rows();i++)
  //      for(int j=0;j<CTI2.cols();j++)
  //        CTI2(i,j)=CTI2(i,j)==-1?-1:(CTI2(i,j)+3-1)%3;
  //  }

  //  assert((CTF.array() == CTF2.array()).all());
  //  assert((CTI.array() == CTI2.array()).all());
  //}
  //std::cout<<igl::matlab_format_index(CF,"CF")<<std::endl;


  SVI.resize(F.size());
  std::vector<bool> marked(F.size());
  VectorXI J = VectorXI::Constant(F.size(),-1);
  SF.resize(F.rows(),F.cols());
  {
    int nv = 0;
    for(int f = 0;f<m;f++)
    {
      for(int i = 0;i<3;i++)
      {
        const int c = CF(f,i);
        if(J(c) == -1)
        {
          J(c) = nv;
          SVI(nv) = F(f,i);
          nv++;
        }
        SF(f,i) = J(c);
      }
    }
    SVI.conservativeResize(nv);
  }

}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedSV,
  typename DerivedSF,
  typename DerivedSVI
  >
IGL_INLINE void igl::split_nonmanifold(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  Eigen::PlainObjectBase <DerivedSV> & SV,
  Eigen::PlainObjectBase <DerivedSF> & SF,
  Eigen::PlainObjectBase <DerivedSVI> & SVI)
{
  igl::split_nonmanifold(F,SF,SVI);
  SV = V(SVI.derived(),Eigen::all);
}

#ifdef IGL_STATIC_LIBRARY
template void igl::split_nonmanifold<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
