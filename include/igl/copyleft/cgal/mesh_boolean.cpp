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
#include "assign.h"
#include "extract_cells.h"
#include "mesh_boolean_type_to_funcs.h"
#include "propagate_winding_numbers.h"
#include "relabel_small_immersed_cells.h"
#include "remesh_self_intersections.h"
#include "string_to_mesh_boolean_type.h"
#include "../../cumsum.h"
#include "../../extract_manifold_patches.h"
#include "../../get_seconds.h"
#include "../../remove_unreferenced.h"
#include "../../resolve_duplicated_faces.h"
#include "../../slice.h"
#include "../../unique_edge_map.h"
#include "../../unique_simplices.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <algorithm>

//#define MESH_BOOLEAN_TIMING
//#define DOUBLE_CHECK_EXACT_OUTPUT
//#define SMALL_CELL_REMOVAL

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA > & VA,
    const Eigen::PlainObjectBase<DerivedFA > & FA,
    const Eigen::PlainObjectBase<DerivedVB > & VB,
    const Eigen::PlainObjectBase<DerivedFB > & FB,
    const MeshBooleanType & type,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J)
{
  std::function<int(const int, const int)> keep;
  std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) > wind_num_op;
  mesh_boolean_type_to_funcs(type,wind_num_op,keep);
  return mesh_boolean(VA,FA,VB,FB,wind_num_op,keep,VC,FC,J);
}
template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const std::string & type_str,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC,
  Eigen::PlainObjectBase<DerivedJ > & J)
{
  return mesh_boolean(
    VA,FA,VB,FB,string_to_mesh_boolean_type(type_str),VC,FC,J);
}

template <
  typename DerivedVA,
  typename DerivedFA,
  typename DerivedVB,
  typename DerivedFB,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVA> & VA,
    const Eigen::PlainObjectBase<DerivedFA> & FA,
    const Eigen::PlainObjectBase<DerivedVB> & VB,
    const Eigen::PlainObjectBase<DerivedFB> & FB,
    const std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) >& wind_num_op,
    const std::function<int(const int, const int)> & keep,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J) 
{
  // Generate combined mesh (VA,FA,VB,FB) -> (V,F)
  Eigen::Matrix<size_t,2,1> sizes(FA.rows(),FB.rows());
  // TODO: This is a precision template **bug** that results in failure to
  // compile. If DerivedVA::Scalar is double and DerivedVB::Scalar is
  // CGAL::Epeck::FT then the following assignment will not compile. This
  // implies that VA must have the trumping precision (and a valid assignment
  // operator from VB's type).
  Eigen::Matrix<typename DerivedVA::Scalar,Eigen::Dynamic,3> VV(VA.rows() + VB.rows(), 3);
  DerivedFC FF(FA.rows() + FB.rows(), 3);
  // Can't use comma initializer
  for(int a = 0;a<VA.rows();a++)
  {
    for(int d = 0;d<3;d++) VV(a,d) = VA(a,d);
  }
  for(int b = 0;b<VB.rows();b++)
  {
    for(int d = 0;d<3;d++) VV(VA.rows()+b,d) = VB(b,d);
  }
  FF.block(0, 0, FA.rows(), 3) = FA;
  FF.block(FA.rows(), 0, FB.rows(), 3) = FB.array() + VA.rows();
  return mesh_boolean(VV,FF,sizes,wind_num_op,keep,VC,FC,J);
}

template <
  typename DerivedV,
  typename DerivedF,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
    const std::vector<DerivedV > & Vlist,
    const std::vector<DerivedF > & Flist,
    const std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) >& wind_num_op,
    const std::function<int(const int, const int)> & keep,
    Eigen::PlainObjectBase<DerivedVC > & VC,
    Eigen::PlainObjectBase<DerivedFC > & FC,
    Eigen::PlainObjectBase<DerivedJ > & J)
{
  assert(Flist.size() == Vlist.size() && "#Vlist and #Flist should match");
  const size_t num_inputs = Vlist.size();
  // Gather sizes
  Eigen::Matrix<size_t,Eigen::Dynamic,1> sizes(num_inputs);
  int numf = 0;
  int numv = 0;
  for(int i = 0;i<num_inputs;i++)
  {
    sizes(i) = Flist[i].rows();
    numf += Flist[i].rows();
    numv += Vlist[i].rows();
  }
  // Combined mesh
  DerivedV VV(numv,3);
  DerivedF FF(numf,3);
  {
    int fk = 0;
    int vk = 0;
    for(int i = 0;i<num_inputs;i++)
    {
      FF.block(fk,0,Flist[i].rows(),3) = Flist[i].array() + vk;
      fk += Flist[i].rows();
      VV.block(vk,0,Vlist[i].rows(),3) = Vlist[i];
      vk += Vlist[i].rows();
    }
  }
  return mesh_boolean(VV,FF,sizes,wind_num_op,keep,VC,FC,J);
}

template <
  typename DerivedVV,
  typename DerivedFF,
  typename Derivedsizes,
  typename DerivedVC,
  typename DerivedFC,
  typename DerivedJ>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
    const Eigen::PlainObjectBase<DerivedVV > & VV,
    const Eigen::PlainObjectBase<DerivedFF > & FF,
    const Eigen::PlainObjectBase<Derivedsizes> & sizes,
    const std::function<int(const Eigen::Matrix<int,1,Eigen::Dynamic>) >& wind_num_op,
    const std::function<int(const int, const int)> & keep,
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
  typedef CGAL::Epeck Kernel;
  typedef Kernel::FT ExactScalar;
  typedef Eigen::Matrix<Scalar,Eigen::Dynamic,3> MatrixX3S;
  typedef Eigen::Matrix<typename DerivedJ::Scalar,Eigen::Dynamic,1> VectorXJ;
  typedef Eigen::Matrix<
    ExactScalar,
    Eigen::Dynamic,
    Eigen::Dynamic,
    DerivedVC::IsRowMajor> MatrixXES;
  MatrixXES V;
  DerivedFC F;
  VectorXJ  CJ;
  {
    Eigen::VectorXi I;
    igl::copyleft::cgal::RemeshSelfIntersectionsParam params;
    params.stitch_all = true;
    MatrixXES Vr;
    DerivedFC Fr;
    Eigen::MatrixXi IF;
    igl::copyleft::cgal::remesh_self_intersections(
        VV, FF, params, Vr, Fr, IF, CJ, I);
    assert(I.size() == Vr.rows());
    // Merge coinciding vertices into non-manifold vertices.
    std::for_each(Fr.data(), Fr.data()+Fr.size(),
          [&I](typename DerivedFC::Scalar& a) { a=I[a]; });
      // Remove unreferenced vertices.
      Eigen::VectorXi UIM;
      igl::remove_unreferenced(Vr, Fr, V, F, UIM);
   }
#ifdef MESH_BOOLEAN_TIMING
  log_time("resolve_self_intersection");
#endif

  // Compute edges of (F) --> (E,uE,EMAP,uE2E)
  Eigen::MatrixXi E, uE;
  Eigen::VectorXi EMAP;
  std::vector<std::vector<size_t> > uE2E;
  igl::unique_edge_map(F, E, uE, EMAP, uE2E);

  // Compute patches (F,EMAP,uE2E) --> (P)
  Eigen::VectorXi P;
  const size_t num_patches = igl::extract_manifold_patches(F, EMAP, uE2E, P);
#ifdef MESH_BOOLEAN_TIMING
  log_time("patch_extraction");
#endif

  // Compute cells (V,F,P,E,uE,EMAP) -> (per_patch_cells)
  Eigen::MatrixXi per_patch_cells;
  const size_t num_cells =
    igl::copyleft::cgal::extract_cells(
        V, F, P, E, uE, uE2E, EMAP, per_patch_cells);
#ifdef MESH_BOOLEAN_TIMING
  log_time("cell_extraction");
#endif

  // Compute winding numbers on each side of each facet.
  const size_t num_faces = F.rows();
  // W(f,:) --> [w1out,w1in,w2out,w2in, ... wnout,wnint] winding numbers above
  // and below each face w.r.t. each input mesh, so that W(f,2*i) is the
  // winding number above face f w.r.t. input i, and W(f,2*i+1) is the winding
  // number below face f w.r.t. input i.
  Eigen::MatrixXi W;
  // labels(f) = i means that face f comes from mesh i
  Eigen::VectorXi labels(num_faces);
  // cumulative sizes
  Derivedsizes cumsizes;
  igl::cumsum(sizes,1,cumsizes);
  const size_t num_inputs = sizes.size();
  std::transform(
    CJ.data(), 
    CJ.data()+CJ.size(), 
    labels.data(),
    // Determine which input mesh birth face i comes from
    [&num_inputs,&cumsizes](int i)->int
    { 
      for(int k = 0;k<num_inputs;k++)
      {
        if(i<cumsizes(k)) return k;
      }
      assert(false && "Birth parent index out of range");
      return -1;
    });
  bool valid = true;
  if (num_faces > 0) 
  {
    valid = valid & 
      igl::copyleft::cgal::propagate_winding_numbers(
          V, F, uE, uE2E, num_patches, P, num_cells, per_patch_cells, labels, W);
  } else 
  {
    W.resize(0, 2*num_inputs);
  }
  assert((size_t)W.rows() == num_faces);
  // If W doesn't have enough columns, pad with zeros
  if (W.cols() <= 2*num_inputs) 
  {
    const int old_ncols = W.cols();
    W.conservativeResize(num_faces,2*num_inputs);
    W.rightCols(2*num_inputs-old_ncols).setConstant(0);
  }
  assert((size_t)W.cols() == 2*num_inputs);
#ifdef MESH_BOOLEAN_TIMING
  log_time("propagate_input_winding_number");
#endif

  // Compute resulting winding number.
  Eigen::MatrixXi Wr(num_faces, 2);
  for (size_t i=0; i<num_faces; i++) 
  {
    // Winding number vectors above and below
    Eigen::RowVectorXi w_out(1,num_inputs), w_in(1,num_inputs);
    for(size_t k =0;k<num_inputs;k++)
    {
      w_out(k) = W(i,2*k+0);
      w_in(k) = W(i,2*k+1);
    }
    Wr(i,0) = wind_num_op(w_out);
    Wr(i,1) = wind_num_op(w_in);
  }
#ifdef MESH_BOOLEAN_TIMING
  log_time("compute_output_winding_number");
#endif

#ifdef SMALL_CELL_REMOVAL
  igl::copyleft::cgal::relabel_small_immersed_cells(
    V, F, num_patches, P, num_cells, per_patch_cells, 1e-3, Wr);
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

    MatrixX3S Vs;
    assign(V,Vs);
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
  typename DerivedVC,
  typename DerivedFC>
IGL_INLINE bool igl::copyleft::cgal::mesh_boolean(
  const Eigen::PlainObjectBase<DerivedVA > & VA,
  const Eigen::PlainObjectBase<DerivedFA > & FA,
  const Eigen::PlainObjectBase<DerivedVB > & VB,
  const Eigen::PlainObjectBase<DerivedFB > & FB,
  const MeshBooleanType & type,
  Eigen::PlainObjectBase<DerivedVC > & VC,
  Eigen::PlainObjectBase<DerivedFC > & FC) 
{
  Eigen::Matrix<typename DerivedFC::Index, Eigen::Dynamic,1> J;
  return igl::copyleft::cgal::mesh_boolean(VA,FA,VB,FB,type,VC,FC,J);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template specialization
// generated by autoexplicit.sh
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3>, Eigen::Matrix<int, 12, 3, 0, 12, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, 12, 3, 0, 12, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
// generated by autoexplicit.sh
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3>, Eigen::Matrix<int, 12, 3, 0, 12, 3>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 4, 0, -1, 4>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 4, 0, -1, 4>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, 12, 3, 0, 12, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 4, 0, -1, 4> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 4, 0, -1, 4> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<long, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<long, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(std::vector<Eigen::Matrix<double, -1, -1, 0, -1, -1>, std::allocator<Eigen::Matrix<double, -1, -1, 0, -1, -1> > > const&, std::vector<Eigen::Matrix<int, -1, -1, 0, -1, -1>, std::allocator<Eigen::Matrix<int, -1, -1, 0, -1, -1> > > const&, std::function<int (Eigen::Matrix<int, 1, -1, 1, 1, -1>)> const&, std::function<int (int, int)> const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3>, Eigen::Matrix<int, 12, 3, 0, 12, 3>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, 8, 3, 0, 8, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, 12, 3, 0, 12, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, std::basic_string<char, std::char_traits<char>, std::allocator<char> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
template bool igl::copyleft::cgal::mesh_boolean<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, igl::MeshBooleanType const&, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&);
#undef IGL_STATIC_LIBRARY
#include "../../remove_unreferenced.cpp"
template void igl::remove_unreferenced<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> >&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> >&);
#include "../../slice.cpp"
template void igl::slice<Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >, Eigen::Matrix<int, -1, 1, 0, -1, 1>, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > >(Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> > const&, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const&, int, Eigen::PlainObjectBase<Eigen::Matrix<CGAL::Lazy_exact_nt<CGAL::Gmpq>, -1, 3, 0, -1, 3> >&);
#endif
