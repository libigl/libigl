// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#include "mesh_boolean.h"
#include <igl/cgal/assign_scalar.h>
#include <igl/per_face_normals.h>
#include <igl/boundary_facets.h>
#include <igl/exterior_edges.h>
#include <igl/cgal/peel_outer_hull_layers.h>
#include <igl/cgal/remesh_self_intersections.h>
#include <igl/remove_unreferenced.h>
#include <igl/mod.h>
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
  typename DerivedJ,
  typename DerivedI>
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
    Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&,
    Eigen::Matrix<typename DerivedI::Scalar, Eigen::Dynamic,1>&)> &
    resolve_fun,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J,
  Eigen::PlainObjectBase<DerivedI > & I)
{
  using namespace Eigen;
  using namespace std;
  using namespace igl;
  using namespace igl::cgal;
  MeshBooleanType eff_type = type;
  // Concatenate A and B into a single mesh
  typedef CGAL::Epeck Kernel;
  typedef Kernel::FT ExactScalar;
  typedef typename DerivedVC::Scalar Scalar;
  typedef typename DerivedFC::Scalar Index;
  typedef Matrix<Scalar,Dynamic,3> MatrixX3S;
  typedef Matrix<ExactScalar,Dynamic,3> MatrixX3ES;
  typedef Matrix<Index,Dynamic,3> MatrixX3I;
  typedef Matrix<Index,Dynamic,2> MatrixX2I;
  typedef Matrix<typename DerivedI::Scalar,Dynamic,1> VectorXI;
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
    VectorXJ & J,
    VectorXI & I)
  {
    MatrixX3ES SV;
    MatrixX3I SF;
    MatrixX2I SIF;
    igl::cgal::RemeshSelfIntersectionsParam params;
    remesh_self_intersections(V,F,params,SV,SF,SIF,J,I);
    for_each(SF.data(),SF.data()+SF.size(),
      [&I](typename MatrixX3I::Scalar & a){a=I(a);});
    {
      Eigen::Matrix<typename MatrixX3S::Index,Dynamic,1> UIM;
      remove_unreferenced(SV,SF,CV,CF,UIM);
      for_each(I.data(),I.data()+I.size(),
        [&UIM](typename VectorXI::Scalar & a){a=UIM(a);});
    }
#ifdef IGL_MESH_BOOLEAN_DEBUG
    cout<<"#F:  "<<F.rows()<<endl;
    cout<<"#CF: "<<CF.rows()<<endl;
#endif
  };

#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"resolve..."<<endl;
#endif
  MatrixX3S CV;
  MatrixX3ES EV;
  MatrixX3I CF;
  VectorXJ CJ;
  Eigen::Matrix<typename DerivedI::Scalar,Dynamic, 1> CI;
  if(resolve_fun)
  {
    resolve_fun(V,F,CV,CF,CJ,CI);
  }else
  {
    libigl_resolve(V,F,EV,CF,CJ,CI);
    CV.resize(EV.rows(), EV.cols());
    // Just use f'ing for loops. What if EV and CV don't use the same ordering?
    for(int i=0;i<EV.rows();i++)
    {
      for(int j=0;j<EV.cols();j++)
      {
        assign_scalar(EV(i,j),CV(i,j));
      }
    }
  }

  if(type == MESH_BOOLEAN_TYPE_RESOLVE)
  {
    FC = CF;
    VC = CV;
    J = CJ;
    I = CI;
    return;
  }

#ifdef IGL_MESH_BOOLEAN_DEBUG
  cout<<"peel..."<<endl;
#endif
  Matrix<bool,Dynamic,1> from_A(CF.rows());
  // peel layers keeping track of odd and even flips
  VectorXi iter;
  Matrix<bool,Dynamic,1> flip;
  peel_outer_hull_layers(EV,CF,iter,flip);
  //Array<bool,Dynamic,1> even = igl::mod(I,2).array()==0;
  const auto even = [&](const Index & f)->bool
  {
    return (iter(f)%2)==0;
  };

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
        if((even(f)&&!flip(f))||(!even(f)&&flip(f)))
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
        if((!even(f) && !flip(f)) || (even(f) && flip(f)))
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
  {
    MatrixXi O;
    boundary_facets(FC,O);
    cout<<"# boundary: "<<O.rows()<<endl;
  }
  cout<<"# exterior: "<<exterior_edges(FC).rows()<<endl;
#endif
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
    VectorXi ucounts = VectorXi::Zero(uG.rows());
    // loop over all faces
    for(Index g = 0;g<gm;g++)
    {
      const int ug = IC(g);
      assert((size_t) ug < uG2G.size());
      uG2G[ug].push_back(g);
      // is uG(g,:) just a rotated version of G(g,:) ?
      const bool consistent =
        (G(g,0) == uG(ug,0) && G(g,1) == uG(ug,1) && G(g,2) == uG(ug,2)) ||
        (G(g,0) == uG(ug,1) && G(g,1) == uG(ug,2) && G(g,2) == uG(ug,0)) ||
        (G(g,0) == uG(ug,2) && G(g,1) == uG(ug,0) && G(g,2) == uG(ug,1));
      counts(ug) += consistent ? 1 : -1;
      ucounts(ug)++;
    }
    MatrixX3I oldG = G;
    // Faces of output vG[i] = j means ith face of output should be jth face in
    // oldG
    vG.clear();
    for(size_t ug = 0;ug < uG2G.size();ug++)
    {
      // if signed occurrences is zero or ±two then keep none
      // else if signed occurrences is ±one then keep just one facet
      switch(abs(counts(ug)))
      {
        case 1:
          assert(uG2G[ug].size() > 0);
          vG.push_back(uG2G[ug][0]);
#ifdef IGL_MESH_BOOLEAN_DEBUG
          if(abs(ucounts(ug)) != 1)
          {
            cout<<"count,ucount of "<<counts(ug)<<","<<ucounts(ug)<<endl;
          }
#endif
          break;
        case 0:
#ifdef IGL_MESH_BOOLEAN_DEBUG
          cout<<"Skipping "<<uG2G[ug].size()<<" facets..."<<endl;
          if(abs(ucounts(ug)) != 0)
          {
            cout<<"count,ucount of "<<counts(ug)<<","<<ucounts(ug)<<endl;
          }
#endif
          break;
        default:
#ifdef IGL_MESH_BOOLEAN_DEBUG
          cout<<"Didn't expect to be here."<<endl;
#endif
          assert(false && "Shouldn't count be -1/0/1 ?");
      }
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
  {
    VectorXi newIM;
    remove_unreferenced(CV,G,VC,FC,newIM);
    I.resize(CI.rows(),CI.cols());
    for(int i = 0;i<CI.rows();i++)
    {
      I(i) = newIM(CI(i));
    }
  }
  //cerr<<"warning not removing unref"<<endl;
  //VC = CV;
  //FC = G;
#ifdef IGL_MESH_BOOLEAN_DEBUG
  {
    MatrixXi O;
    boundary_facets(FC,O);
    cout<<"# boundary: "<<O.rows()<<endl;
  }
  cout<<"# exterior: "<<exterior_edges(FC).rows()<<endl;
#endif
}


template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ,
  typename DerivedI>
IGL_INLINE void igl::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J,
  Eigen::PlainObjectBase<DerivedI > & I)
{
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedJ::Scalar, Eigen::Dynamic,1>&,
          Eigen::Matrix<typename DerivedI::Scalar, Eigen::Dynamic,1>&) >
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J,I);
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
  typedef Eigen::Matrix<typename DerivedVC::Index, Eigen::Dynamic,1> VectorXI;
  VectorXI I;
  const std::function<void(
    const Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
    const Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedVC::Scalar,Eigen::Dynamic,3> &,
          Eigen::Matrix<typename DerivedFC::Scalar, Eigen::Dynamic,3>&,
          Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1>&,
          Eigen::Matrix<typename VectorXI::Scalar, Eigen::Dynamic,1>&)>
    empty_fun;
  return mesh_boolean(VA,FA,VB,FB,type,empty_fun,VC,FC,J,I);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template void igl::boolean::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
// This is a hack to discuss. I'm not sure why this _doesn't_ create
// duplicate symbols.
#include <igl/remove_unreferenced.cpp>
template void igl::remove_unreferenced<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#include <igl/cgal/peel_outer_hull_layers.cpp>
template unsigned long
igl::cgal::peel_outer_hull_layers<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#include <igl/cgal/outer_hull.cpp>
template void igl::cgal::outer_hull<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, Eigen::Matrix<bool, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&, Eigen::PlainObjectBase<Eigen::Matrix<bool, -1, 1, 0, -1, 1> >&);
#include <igl/slice.cpp>
#include <igl/barycenter.cpp>
template void igl::barycenter<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&);
#include <igl/mod.cpp>
template Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > igl::mod<Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int);
#include <igl/outer_element.cpp>
template void igl::outer_edge<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, 3, 0, -1, 3>, Eigen::Matrix<long, -1, 1, 0, -1, 1>, long, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> > const&, long&, long&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
#include <igl/colon.cpp>
template void igl::colon<int, long, long>(int, long, Eigen::Matrix<long, -1, 1, 0, -1, 1>&);


#endif
