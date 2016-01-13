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
#include "../../get_seconds.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <algorithm>

//#define MESH_BOOLEAN_TIMING
//#define DOUBLE_CHECK_EXACT_OUTPUT

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
IGL_INLINE bool igl::copyleft::boolean::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA> & VA,
    const Eigen::PlainObjectBase<DerivedFA> & FA,
    const Eigen::PlainObjectBase<DerivedVB> & VB,
    const Eigen::PlainObjectBase<DerivedFB> & FB,
    const WindingNumberOp& wind_num_op,
    const KeepFunc& keep,
    const ResolveFunc& resolve_fun,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J) 
{

#ifdef MESH_BOOLEAN_TIMING
  const auto & tictoc = []() -> double
  {
    static double t_start = igl::get_seconds();
    double diff = igl::get_seconds()-t_start;
    t_start += diff;
    return diff;
  };
  const auto log_time = [&](const std::string& label) -> void {
    std::cout << "mesh_boolean." << label << ": "
      << tictoc() << std::endl;
  };
  tictoc();
#endif

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
      resolve_fun(VV, FF, V, F, CJ);
  }
#ifdef MESH_BOOLEAN_TIMING
  log_time("resolve_self_intersection");
#endif

  // Compute winding numbers on each side of each facet.
  const size_t num_faces = F.rows();
  Eigen::MatrixXi W;
  Eigen::VectorXi labels(num_faces);
  std::transform(CJ.data(), CJ.data()+CJ.size(), labels.data(),
      [&](int i) { return i<FA.rows() ? 0:1; });
  bool valid = true;
  if (num_faces > 0) 
  {
    valid = valid & 
      igl::copyleft::cgal::propagate_winding_numbers(V, F, labels, W);
  } else 
  {
    W.resize(0, 4);
  }
  assert((size_t)W.rows() == num_faces);
  if (W.cols() == 2) 
  {
    assert(FB.rows() == 0);
    Eigen::MatrixXi W_tmp(num_faces, 4);
    W_tmp << W, Eigen::MatrixXi::Zero(num_faces, 2);
    W = W_tmp;
  } else {
    assert(W.cols() == 4);
  }
#ifdef MESH_BOOLEAN_TIMING
  log_time("propagate_input_winding_number");
#endif

  // Compute resulting winding number.
  Eigen::MatrixXi Wr(num_faces, 2);
  for (size_t i=0; i<num_faces; i++) 
  {
    Eigen::MatrixXi w_out(1,2), w_in(1,2);
    w_out << W(i,0), W(i,2);
    w_in  << W(i,1), W(i,3);
    Wr(i,0) = wind_num_op(w_out);
    Wr(i,1) = wind_num_op(w_in);
  }
#ifdef MESH_BOOLEAN_TIMING
  log_time("compute_output_winding_number");
#endif

  // Extract boundary separating inside from outside.
  auto index_to_signed_index = [&](size_t i, bool ori) -> int
  {
    return (i+1)*(ori?1:-1);
  };
  //auto signed_index_to_index = [&](int i) -> size_t {
  //    return abs(i) - 1;
  //};
  std::vector<int> selected;
  for(size_t i=0; i<num_faces; i++) 
  {
    auto should_keep = keep(Wr(i,0), Wr(i,1));
    if (should_keep > 0) 
    {
      selected.push_back(index_to_signed_index(i, true));
    } else if (should_keep < 0) 
    {
      selected.push_back(index_to_signed_index(i, false));
    }
  }

  const size_t num_selected = selected.size();
  DerivedFC kept_faces(num_selected, 3);
  DerivedJ  kept_face_indices(num_selected, 1);
  for (size_t i=0; i<num_selected; i++) 
  {
    size_t idx = abs(selected[i]) - 1;
    if (selected[i] > 0) 
    {
      kept_faces.row(i) = F.row(idx);
    } else 
    {
      kept_faces.row(i) = F.row(idx).reverse();
    }
    kept_face_indices(i, 0) = CJ[idx];
  }
#ifdef MESH_BOOLEAN_TIMING
  log_time("extract_output");
#endif

  // Finally, remove duplicated faces and unreferenced vertices.
  {
    DerivedFC G;
    DerivedJ JJ;
    igl::resolve_duplicated_faces(kept_faces, G, JJ);
    igl::slice(kept_face_indices, JJ, 1, J);

#ifdef DOUBLE_CHECK_EXACT_OUTPUT
    {
      // Sanity check on exact output.
      igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
      params.detect_only = true;
      params.first_only = true;
      MatrixXES dummy_VV;
      DerivedFC dummy_FF, dummy_IF;
      Eigen::VectorXi dummy_J, dummy_IM;
      igl::copyleft::cgal::SelfIntersectMesh<
        Kernel,
        MatrixXES, DerivedFC,
        MatrixXES, DerivedFC,
        DerivedFC,
        Eigen::VectorXi,
        Eigen::VectorXi
      > checker(V, G, params,
          dummy_VV, dummy_FF, dummy_IF, dummy_J, dummy_IM);
      if (checker.count != 0) 
      {
        throw "Self-intersection not fully resolved.";
      }
    }
#endif

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
#ifdef MESH_BOOLEAN_TIMING
  log_time("clean_up");
#endif
  return valid;
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
IGL_INLINE bool igl::copyleft::boolean::mesh_boolean(
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
  switch (type) 
  {
    case MESH_BOOLEAN_TYPE_UNION:
      return igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryUnion(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
    case MESH_BOOLEAN_TYPE_INTERSECT:
      return igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryIntersect(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
    case MESH_BOOLEAN_TYPE_MINUS:
      return igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryMinus(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
    case MESH_BOOLEAN_TYPE_XOR:
      return igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryXor(),
          igl::copyleft::boolean::KeepInside(), resolve_func, VC, FC, J);
    case MESH_BOOLEAN_TYPE_RESOLVE:
      //op = binary_resolve();
      return igl::copyleft::boolean::mesh_boolean(
          VA, FA, VB, FB, igl::copyleft::boolean::BinaryResolve(),
          igl::copyleft::boolean::KeepAll(), resolve_func, VC, FC, J);
    default:
      assert(false && "Unsupported boolean type.");
      return false;
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
IGL_INLINE bool igl::copyleft::boolean::mesh_boolean(
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
IGL_INLINE bool igl::copyleft::boolean::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC) 
{
  Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1> J;
  return igl::copyleft::boolean::mesh_boolean(VA,FA,VB,FB,type,VC,FC,J);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
template bool igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::boolean::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::copyleft::boolean::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#undef IGL_STATIC_LIBRARY
#include "../../remove_unreferenced.cpp"
template void igl::remove_unreferenced<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#include "../../slice.cpp"
template void igl::slice<Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&);
#endif
