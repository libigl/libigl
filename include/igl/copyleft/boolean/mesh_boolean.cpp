// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//                    Qingnan Zhou <qnzhou@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
//
#include "mesh_boolean.h"
#include "BinaryWindingNumberOperations.h"
#include "../cgal/assign_scalar.h"
#include "../cgal/propagate_winding_numbers.h"
#include "../cgal/remesh_self_intersections.h"
#include "../../remove_unreferenced.h"
#include "../../unique_simplices.h"
#include "../../slice.h"
#include "../../resolve_duplicated_faces.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <algorithm>

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename WindingNumberOp,
  typename KeepFunc,
  typename ResolveFunc,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::copyleft::boolean::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA> & VA,
    const Eigen::PlainObjectBase<DerivedFA> & FA,
    const Eigen::PlainObjectBase<DerivedVB> & VB,
    const Eigen::PlainObjectBase<DerivedFB> & FB,
    const WindingNumberOp& wind_num_op,
    const KeepFunc& keep,
    const ResolveFunc& resolve_fun,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J) {

  typedef typename DerivedVC::Scalar Scalar;
  //typedef typename DerivedFC::Scalar Index;
  typedef CGAL::Epeck Kernel;
  typedef Kernel::FT ExactScalar;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,3> MatrixX3S;
  //typedef Eigen::Matrix<Index,Eigen::Dynamic,Eigen::Dynamic> MatrixXI;
  typedef Eigen::Matrix<typename DerivedJ::Scalar,Eigen::Dynamic,1> VectorXJ;

  // Generate combined mesh.
  typedef Eigen::Matrix<
    ExactScalar,
    Eigen::Dynamic,
    Eigen::Dynamic,
    DerivedVC::IsRowMajor> MatrixXES;
  MatrixXES V;
  DerivedFC F;
  VectorXJ  CJ;
  {
      DerivedVA VV(VA.rows() + VB.rows(), 3);
      DerivedFC FF(FA.rows() + FB.rows(), 3);
      VV << VA, VB;
      FF << FA, FB.array() + VA.rows();
      //// Handle annoying empty cases
      //if(VA.size()>0)
      //{
      //  VV<<VA;
      //}
      //if(VB.size()>0)
      //{
      //  VV<<VB;
      //}
      //if(FA.size()>0)
      //{
      //  FF<<FA;
      //}
      //if(FB.size()>0)
      //{
      //  FF<<FB.array()+VA.rows();
      //}
      resolve_fun(VV, FF, V, F, CJ);
  }

  // Compute winding numbers on each side of each facet.
  const size_t num_faces = F.rows();
  Eigen::MatrixXi W;
  Eigen::VectorXi labels(num_faces);
  std::transform(CJ.data(), CJ.data()+CJ.size(), labels.data(),
      [&](int i) { return i<FA.rows() ? 0:1; });
  igl::copyleft::cgal::propagate_winding_numbers(V, F, labels, W);
  assert((size_t)W.rows() == num_faces);
  if (W.cols() == 2) {
    assert(FB.rows() == 0);
    Eigen::MatrixXi W_tmp(num_faces, 4);
    W_tmp << W, Eigen::MatrixXi::Zero(num_faces, 2);
    W = W_tmp;
  } else {
    assert(W.cols() == 4);
  }

  // Compute resulting winding number.
  Eigen::MatrixXi Wr(num_faces, 2);
  for (size_t i=0; i<num_faces; i++) {
    Eigen::MatrixXi w_out(1,2), w_in(1,2);
    w_out << W(i,0), W(i,2);
    w_in  << W(i,1), W(i,3);
    Wr(i,0) = wind_num_op(w_out);
    Wr(i,1) = wind_num_op(w_in);
  }

  // Extract boundary separating inside from outside.
  auto index_to_signed_index = [&](size_t i, bool ori) -> int{
    return (i+1)*(ori?1:-1);
  };
  //auto signed_index_to_index = [&](int i) -> size_t {
  //    return abs(i) - 1;
  //};
  std::vector<int> selected;
  for(size_t i=0; i<num_faces; i++) {
    auto should_keep = keep(Wr(i,0), Wr(i,1));
    if (should_keep > 0) {
      selected.push_back(index_to_signed_index(i, true));
    } else if (should_keep < 0) {
      selected.push_back(index_to_signed_index(i, false));
    }
  }

  const size_t num_selected = selected.size();
  DerivedFC kept_faces(num_selected, 3);
  DerivedJ  kept_face_indices(num_selected, 1);
  for (size_t i=0; i<num_selected; i++) {
    size_t idx = abs(selected[i]) - 1;
    if (selected[i] > 0) {
      kept_faces.row(i) = F.row(idx);
    } else {
      kept_faces.row(i) = F.row(idx).reverse();
    }
    kept_face_indices(i, 0) = CJ[idx];
  }

  // Finally, remove duplicated faces and unreferenced vertices.
  {
    DerivedFC G;
    DerivedJ JJ;
    igl::resolve_duplicated_faces(kept_faces, G, JJ);
    igl::slice(kept_face_indices, JJ, 1, J);

    MatrixX3S Vs(V.rows(), V.cols());
    for (size_t i=0; i<(size_t)V.rows(); i++)
    {
      for (size_t j=0; j<(size_t)V.cols(); j++)
      {
        igl::copyleft::cgal::assign_scalar(V(i,j), Vs(i,j));
      }
    }
    Eigen::VectorXi newIM;
    igl::remove_unreferenced(Vs,G,VC,FC,newIM);
  }
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename ResolveFunc,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::copyleft::boolean::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA > & VA,
    const Eigen::PlainObjectBase<DerivedFA > & FA,
    const Eigen::PlainObjectBase<DerivedVB > & VB,
    const Eigen::PlainObjectBase<DerivedFB > & FB,
    const MeshBooleanType & type,
    const ResolveFunc& resolve_func,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J)
{
  //typedef CGAL::Epeck Kernel;
  //typedef Kernel::FT ExactScalar;
  //typedef Eigen::Matrix<
  //  ExactScalar,
  //  Eigen::Dynamic,
  //  Eigen::Dynamic,
  //  DerivedVC::IsRowMajor> MatrixXES;

  switch (type) {
    case MESH_BOOLEAN_TYPE_UNION:
      igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryUnion(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
      break;
    case MESH_BOOLEAN_TYPE_INTERSECT:
      igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryIntersect(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
      break;
    case MESH_BOOLEAN_TYPE_MINUS:
      igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryMinus(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
      break;
    case MESH_BOOLEAN_TYPE_XOR:
      igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryXor(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
      break;
    case MESH_BOOLEAN_TYPE_RESOLVE:
      //op = binary_resolve();
      igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryResolve(),
          igl::copyleft::boolean::KeepAll(), resolve_func, VC, FC, J);
      break;
    default:
      throw std::runtime_error("Unsupported boolean type.");
  }
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE void igl::copyleft::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  typedef CGAL::Epeck Kernel;
  typedef Kernel::FT ExactScalar;
  typedef Eigen::Matrix<
    ExactScalar,
    Eigen::Dynamic,
    Eigen::Dynamic,
    DerivedVC::IsRowMajor> MatrixXES;

  std::function<void(
      const Eigen::PlainObjectBase<DerivedVA>&,
      const Eigen::PlainObjectBase<DerivedFA>&,
      Eigen::PlainObjectBase<MatrixXES>&,
      Eigen::PlainObjectBase<DerivedFC>&,
      Eigen::PlainObjectBase<DerivedJ>&)> resolve_func =
    [](const Eigen::PlainObjectBase<DerivedVA>& V,
        const Eigen::PlainObjectBase<DerivedFA>& F,
        Eigen::PlainObjectBase<MatrixXES>& Vo,
        Eigen::PlainObjectBase<DerivedFC>& Fo,
        Eigen::PlainObjectBase<DerivedJ>& J) {
      Eigen::VectorXi I;
      igl::copyleft::cgal::RemeshSelfIntersectionsParam params;

      MatrixXES Vr;
      DerivedFC Fr;
      Eigen::MatrixXi IF;
      igl::copyleft::cgal::remesh_self_intersections(
          V, F, params, Vr, Fr, IF, J, I);
      assert(I.size() == Vr.rows());

      // Merge coinciding vertices into non-manifold vertices.
      std::for_each(Fr.data(), Fr.data()+Fr.size(),
          [&I](typename DerivedFC::Scalar& a) { a=I[a]; });

      // Remove unreferenced vertices.
      Eigen::VectorXi UIM;
      igl::remove_unreferenced(Vr, Fr, Vo, Fo, UIM);
    };

  return mesh_boolean(VA,FA,VB,FB,type,resolve_func,VC,FC,J);
}

template <
typename DerivedVA,
typename DerivedFA,
typename DerivedVB,
typename DerivedFB,
typename DerivedVC,
typename DerivedFC>
IGL_INLINE void igl::copyleft::boolean::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA > & VA,
    const Eigen::PlainObjectBase<DerivedFA > & FA,
    const Eigen::PlainObjectBase<DerivedVB > & VB,
    const Eigen::PlainObjectBase<DerivedFB > & FB,
    const MeshBooleanType & type,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC) {
  Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1> J;
  return igl::copyleft::boolean::mesh_boolean(VA,FA,VB,FB,type,VC,FC,J);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template void igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template void igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template void igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#undef IGL_STATIC_LIBRARY
#include "../../remove_unreferenced.cpp"
template void igl::remove_unreferenced<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#include "../../slice.cpp"
template void igl::slice<Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&);
#endif
