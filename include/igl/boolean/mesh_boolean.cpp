// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mesh_boolean.h"
#include <igl/per_face_normals.h>
#include <igl/cgal/peel_outer_hull_layers.h>
#include <igl/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>
#include <igl/unique_simplices.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <iostream>

//#define IGL_MESH_BOOLEAN_DEBUG

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&)>
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J);
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC>
IGL_INLINE void igl::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC)
{
  Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1> J;
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1>&)>
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J);
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3>&,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&)>
    & resolve_fun,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  using namespace igl::cgal;
  MeshBooleanType eff_type = type;
  // Concatenate A and B into a single mesh
  typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
  typedef Kernel::FT ExactScalar;
  typedef typename DerivedVC::Scalar Scalar;
  typedef typename DerivedFC::Scalar Index;
  typedef Matrix<Scalar,Dynamic,3> MatrixX3S;
  typedef Matrix<ExactScalar,Dynamic,3> MatrixX3ES;
  typedef Matrix<Index,Dynamic,3> MatrixX3I;
  typedef Matrix<Index,Dynamic,2> MatrixX2I;
  typedef Matrix<Index,Dynamic,1> VectorXI;
  typedef Matrix<typename DerivedJ::Scalar,Dynamic,1> VectorXJ;
#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"mesh boolean..."<<endl;
#endif
  MatrixX3S V(VA.rows()+VB.rows(),3);
  MatrixX3I F(FA.rows()+FB.rows(),3);
  V.block(0,0,VA.rows(),VA.cols()) = VA;
  V.block(VA.rows(),0,VB.rows(),VB.cols()) = VB;
#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"prepare selfintersect input..."<<endl;
#endif
  switch(type)
  {
    // Minus is implemented by flipping B and computing union
    case MESH_BOOLEAN_TYPE_MINUS:
      F.block(0,0,FA.rows(),FA.cols()) = FA.rowwise().reverse();
      F.block(FA.rows(),0,FB.rows(),FB.cols()) = FB.array()+VA.rows();
      //F.block(0,0,FA.rows(),3) = FA;
      //F.block(FA.rows(),0,FB.rows(),3) =
      //  FB.rowwise().reverse().array()+VA.rows();
      eff_type = MESH_BOOLEAN_TYPE_INTERSECT;
      break;
    default:
      F.block(0,0,FA.rows(),FA.cols()) = FA;
      F.block(FA.rows(),0,FB.rows(),FB.cols()) = FB.array()+VA.rows();
      break;
  }

  // Resolve intersections (assumes A and B are solid)
  const auto & libigl_resolve = [](
    const MatrixX3S & V,
    const MatrixX3I & F,
    MatrixX3ES & CV,
    MatrixX3I & CF,
    VectorXJ & J)
  {
    MatrixX3ES SV;
    MatrixX3I SF;
    MatrixX2I SIF;
    VectorXI SIM,UIM;
    igl::cgal::RemeshSelfIntersectionsParam params;
    remesh_self_intersections(V,F,params,SV,SF,SIF,J,SIM);
    for_each(SF.data(),SF.data()+SF.size(),[&SIM](int & a){a=SIM(a);});
    {
      remove_unreferenced(SV,SF,CV,CF,UIM);
    }
  };

#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"resolve..."<<endl;
#endif
  MatrixX3S CV;
  MatrixX3ES EV;
  MatrixX3I CF;
  VectorXJ CJ;
  if(resolve_fun)
  {
    resolve_fun(V,F,CV,CF,CJ);
  }else
  {
    libigl_resolve(V,F,EV,CF,CJ);
    CV.resize(EV.rows(), EV.cols());
    std::transform(EV.data(), EV.data() + EV.rows()*EV.cols(),
            CV.data(), [&](ExactScalar val) {
            return CGAL::to_double(val);
            });
  }

  if(type == MESH_BOOLEAN_TYPE_RESOLVE)
  {
    FC = CF;
    VC = CV;
    J = CJ;
    return;
  }
  MatrixX3S N,CN;
  per_face_normals_stable(V,F,N);
  CN.resize(CF.rows(),3);
  for(size_t f = 0;f<(size_t)CN.rows();f++)
  {
    CN.row(f) = N.row(CJ(f));
  }

#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"peel..."<<endl;
#endif
  Matrix<bool,Dynamic,1> from_A(CF.rows());
  // peel layers keeping track of odd and even flips
  Matrix<bool,Dynamic,1> odd;
  Matrix<bool,Dynamic,1> flip;
  peel_outer_hull_layers(EV,CF,CN,odd,flip);

#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"categorize..."<<endl;
#endif
  const Index m = CF.rows();
  // Faces of output vG[i] = j means ith face of output should be jth face in F
  std::vector<Index> vG;
  // Whether faces of output should be flipped, Gflip[i] = true means ith face
  // of output should be F.row(vG[i]).reverse() rather than F.row(vG[i])
  std::vector<bool> Gflip;
  for(Index f = 0;f<m;f++)
  {
    switch(eff_type)
    {
      case MESH_BOOLEAN_TYPE_XOR:
      case MESH_BOOLEAN_TYPE_UNION:
        if((odd(f)&&!flip(f))||(!odd(f)&&flip(f)))
        {
          vG.push_back(f);
          Gflip.push_back(false);
        }else if(eff_type == MESH_BOOLEAN_TYPE_XOR)
        {
          vG.push_back(f);
          Gflip.push_back(true);
        }
        break;
      case MESH_BOOLEAN_TYPE_INTERSECT:
        if((!odd(f) && !flip(f)) || (odd(f) && flip(f)))
        {
          vG.push_back(f);
          Gflip.push_back(type == MESH_BOOLEAN_TYPE_MINUS);
        }
        break;
      default:
        assert(false && "Unknown type");
        return;
    }
  }
  const Index gm = vG.size();
  MatrixX3I G(gm,3);
  VectorXi GJ(gm,1);
  for(Index g = 0;g<gm;g++)
  {
    G.row(g) = Gflip[g] ? CF.row(vG[g]).reverse().eval() : CF.row(vG[g]);
    GJ(g) = CJ(vG[g]);
  }
#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"clean..."<<endl;
#endif
  // Deal with duplicate faces
  {
    VectorXi IA,IC;
    MatrixX3I uG;
    unique_simplices(G,uG,IA,IC);
    assert(IA.rows() == uG.rows());
    // faces ontop of each unique face
    vector<vector<Index> > uG2G(uG.rows());
    // signed counts
    VectorXi counts = VectorXi::Zero(uG.rows());
    // loop over all faces
    for(Index g = 0;g<gm;g++)
    {
      const int ug = IC(g);
      assert(ug < uG2G.size());
      uG2G[ug].push_back(g);
      // is uG(g,:) just a rotated version of G(g,:) ?
      const bool consistent =
        (G(g,0) == uG(ug,0) && G(g,1) == uG(ug,1) && G(g,2) == uG(ug,2)) ||
        (G(g,0) == uG(ug,1) && G(g,1) == uG(ug,2) && G(g,2) == uG(ug,0)) ||
        (G(g,0) == uG(ug,2) && G(g,1) == uG(ug,0) && G(g,2) == uG(ug,1));
      counts(ug) += consistent ? 1 : -1;
    }
    MatrixX3I oldG = G;
    // Faces of output vG[i] = j means ith face of output should be jth face in
    // oldG
    vG.clear();
    for(size_t ug = 0;ug < uG2G.size();ug++)
    {
      // if signed occurrences is zero or ±two then keep none
      // else if signed occurrences is ±one then keep just one facet
      if(abs(counts(ug)) == 1)
      {
        assert(uG2G.size() > 0);
        vG.push_back(uG2G[ug][0]);
      }
#ifdef IGL_MESH_BOOLEAN_DEBUG
      else
      {
        cout<<"Skipping "<<uG2G.size()<<" facets..."<<endl;
      }
#endif
    }
    G.resize(vG.size(),3);
    J.resize(vG.size());
    for(size_t g = 0;g<vG.size();g++)
    {
      G.row(g) = oldG.row(vG[g]);
      J(g) = GJ(vG[g]);
    }
  }
  // remove unreferenced vertices
  VectorXi newIM;
  remove_unreferenced(CV,G,VC,FC,newIM);
  //cerr<<"warning not removing unref"<<endl;
  //VC = CV;
  //FC = G;
}

#ifdef IGL_STATIC_LIBRARY
// This is a hack to discuss. I'm not sure why this _doesn't_ create
// duplicate symbols.
#include <igl/remove_unreferenced.cpp>
template void igl::remove_unreferenced<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#include <igl/cgal/peel_outer_hull_layers.cpp>
template unsigned long igl::cgal::peel_outer_hull_layers<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<bool, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#include <igl/cgal/outer_hull.cpp>
template void igl::cgal::outer_hull<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<double, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
// Explicit template specialization
template void igl::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#endif
