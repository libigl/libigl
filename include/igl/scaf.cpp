// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "scaf.h"

#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <igl/ARAPEnergyType.h>
#include <igl/PI.h>
#include <igl/Timer.h>
#include <igl/Timer.h>
#include <igl/arap.h>
#include <igl/arap_linear_block.h>
#include <igl/boundary_facets.h>
#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/cat.h>
#include <igl/colon.h>
#include <igl/cotmatrix.h>
#include <igl/covariance_scatter_matrix.h>
#include <igl/doublearea.h>
#include <igl/edge_flaps.h>
#include <igl/edge_lengths.h>
#include <igl/flip_avoiding_line_search.h>
#include <igl/flipped_triangles.h>
#include <igl/grad.h>
#include <igl/harmonic.h>
#include <igl/is_symmetric.h>
#include <igl/local_basis.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/massmatrix.h>
#include <igl/per_face_normals.h>
#include <igl/polar_svd.h>
#include <igl/readOBJ.h>
#include <igl/repdiag.h>
#include <igl/serialize.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/slim.h>
#include <igl/triangle/triangulate.h>
#include <igl/unique.h>
#include <igl/vector_area_matrix.h>
#include <igl/volume.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <map>
#include <algorithm>
#include <set>
#include <vector>

namespace igl
{
namespace scaf
{
void update_scaffold(igl::SCAFData &s)
{
  s.mv_num = s.m_V.rows();
  s.mf_num = s.m_T.rows();

  s.v_num = s.w_uv.rows();
  s.sf_num = s.s_T.rows();

  s.sv_num = s.v_num - s.mv_num;
  s.f_num = s.sf_num + s.mf_num;

  s.s_M = Eigen::VectorXd::Constant(s.sf_num, s.scaffold_factor);
}

void add_soft_constraints(igl::SCAFData &s, const Eigen::VectorXi &b,
                          const Eigen::MatrixXd &bc)
{
  assert(b.rows() == bc.rows() && "Constraint input incompatible");
  for (int i = 0; i < b.rows(); i++)
    s.soft_cons[b(i)] = bc.row(i);
}

void add_soft_constraints(igl::SCAFData &s, int b, const Eigen::RowVectorXd &bc)
{
  s.soft_cons[b] = bc;
}

void automatic_expand_frame(igl::SCAFData &s, double min2, double max3)
{
  // right top
  // left down
  using namespace Eigen;
  MatrixXd m_uv = s.w_uv.topRows(s.mv_num);
  MatrixXd frame(2, s.dim), bbox(2, s.dim);
  frame << s.w_uv.colwise().maxCoeff(), s.w_uv.colwise().minCoeff();
  bbox << m_uv.colwise().maxCoeff(), m_uv.colwise().minCoeff();
  RowVector2d center = bbox.colwise().mean();
  /*
        bbox.row(0) -= center;
        bbox.row(1) -= center;
        frame.row(0) -= center;
        frame.row(1) -= center;
      */
  struct line_func
  {
    double a, b;

    double operator()(double y)
    {
      return a * y + b;
    };
  };

  auto linear_stretch = [](double s0,
                           double t0,
                           double s1,
                           double t1) { // source0, target0, source1, target1
    Matrix2d S;
    S << s0, 1, s1, 1;
    Vector2d t;
    t << t0, t1;
    Vector2d coef = S.colPivHouseholderQr().solve(t);
    return line_func{coef[0], coef[1]};
  };

  double new_frame;
  double center_coord;
  for (auto d = 0; d < s.dim; d++)
  {
    center_coord = center(d);

    if (frame(0, d) - center_coord < min2 * (bbox(0, d) - center_coord))
    {
      new_frame = max3 * (bbox(0, d) - center_coord) + center_coord;
      auto expand = linear_stretch(bbox(0, d), bbox(0, d),
                                   frame(0, d), new_frame);
      for (auto v = s.mv_num; v < s.v_num; v++)
      {
        if (s.w_uv(v, d) > bbox(0, d))
          s.w_uv(v, d) = expand(s.w_uv(v, d));
      }
    }

    if (frame(1, d) - center_coord > min2 * (bbox(1, d) - center_coord))
    {
      new_frame = max3 * (bbox(1, d) - center_coord) + center_coord;
      auto expand = linear_stretch(bbox(1, d), bbox(1, d),
                                   frame(1, d), new_frame);
      for (auto v = s.mv_num; v < s.v_num; v++)
      {
        if (s.w_uv(v, d) < bbox(1, d))
          s.w_uv(v, d) = expand(s.w_uv(v, d));
      }
    }
  }
}

// functions from ReweightedARAP, a static function
void compute_surface_gradient_matrix(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &F1,
    const Eigen::MatrixXd &F2,
    Eigen::SparseMatrix<double> &D1,
    Eigen::SparseMatrix<double> &D2)
{
  Eigen::SparseMatrix<double> G;
  igl::grad(V, F, G);
  Eigen::SparseMatrix<double> Dx = G.block(0, 0, F.rows(), V.rows());
  Eigen::SparseMatrix<double> Dy = G.block(F.rows(), 0, F.rows(), V.rows());
  Eigen::SparseMatrix<double> Dz = G.block(2 * F.rows(), 0, F.rows(), V.rows());

  D1 = F1.col(0).asDiagonal() * Dx + F1.col(1).asDiagonal() * Dy +
       F1.col(2).asDiagonal() * Dz;
  D2 = F2.col(0).asDiagonal() * Dx + F2.col(1).asDiagonal() * Dy +
       F2.col(2).asDiagonal() * Dz;
}

template <typename DerivedV, typename DerivedF>
inline void adjusted_grad(const Eigen::MatrixBase<DerivedV> &V,
                          const Eigen::MatrixBase<DerivedF> &F,
                          Eigen::SparseMatrix<typename DerivedV::Scalar> &G,
                          double eps)
{
  Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 3>
      eperp21(F.rows(), 3), eperp13(F.rows(), 3);
  int fixed = 0;
  for (int i = 0; i < F.rows(); ++i)
  {
    // renaming indices of vertices of triangles for convenience
    int i1 = F(i, 0);
    int i2 = F(i, 1);
    int i3 = F(i, 2);

    // #F x 3 matrices of triangle edge vectors, named after opposite vertices
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v32 = V.row(i3) - V.row(i2);
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v13 = V.row(i1) - V.row(i3);
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> v21 = V.row(i2) - V.row(i1);
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> n = v32.cross(v13);
    // area of parallelogram is twice area of triangle
    // area of parallelogram is || v1 x v2 ||
    // This does correct l2 norm of rows, so that it contains #F list of twice
    // triangle areas
    double dblA = std::sqrt(n.dot(n));
    Eigen::Matrix<typename DerivedV::Scalar, 1, 3> u;
    if (dblA > eps)
    {
      // now normalize normals to get unit normals
      u = n / dblA;
    }
    else
    {
      // Abstract equilateral triangle v1=(0,0), v2=(h,0), v3=(h/2, (sqrt(3)/2)*h)
      fixed++;
      // get h (by the area of the triangle)
      dblA = eps;
      double h = sqrt((dblA) / sin(
                                   M_PI / 3.0)); // (h^2*sin(60))/2. = Area => h = sqrt(2*Area/sin_60)

      Eigen::Vector3d v1, v2, v3;
      v1 << 0, 0, 0;
      v2 << h, 0, 0;
      v3 << h / 2., (sqrt(3) / 2.) * h, 0;

      // now fix v32,v13,v21 and the normal
      v32 = v3 - v2;
      v13 = v1 - v3;
      v21 = v2 - v1;
      n = v32.cross(v13);
    }

    // rotate each vector 90 degrees around normal
    double norm21 = std::sqrt(v21.dot(v21));
    double norm13 = std::sqrt(v13.dot(v13));
    eperp21.row(i) = u.cross(v21);
    eperp21.row(i) =
        eperp21.row(i) / std::sqrt(eperp21.row(i).dot(eperp21.row(i)));
    eperp21.row(i) *= norm21 / dblA;
    eperp13.row(i) = u.cross(v13);
    eperp13.row(i) =
        eperp13.row(i) / std::sqrt(eperp13.row(i).dot(eperp13.row(i)));
    eperp13.row(i) *= norm13 / dblA;
  }

  std::vector<int> rs;
  rs.reserve(F.rows() * 4 * 3);
  std::vector<int> cs;
  cs.reserve(F.rows() * 4 * 3);
  std::vector<double> vs;
  vs.reserve(F.rows() * 4 * 3);

  // row indices
  for (int r = 0; r < 3; r++)
  {
    for (int j = 0; j < 4; j++)
    {
      for (int i = r * F.rows(); i < (r + 1) * F.rows(); i++)
        rs.push_back(i);
    }
  }

  // column indices
  for (int r = 0; r < 3; r++)
  {
    for (int i = 0; i < F.rows(); i++)
      cs.push_back(F(i, 1));
    for (int i = 0; i < F.rows(); i++)
      cs.push_back(F(i, 0));
    for (int i = 0; i < F.rows(); i++)
      cs.push_back(F(i, 2));
    for (int i = 0; i < F.rows(); i++)
      cs.push_back(F(i, 0));
  }

  // values
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp13(i, 0));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp13(i, 0));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp21(i, 0));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp21(i, 0));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp13(i, 1));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp13(i, 1));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp21(i, 1));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp21(i, 1));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp13(i, 2));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp13(i, 2));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(eperp21(i, 2));
  for (int i = 0; i < F.rows(); i++)
    vs.push_back(-eperp21(i, 2));

  // create sparse gradient operator matrix
  G.resize(3 * F.rows(), V.rows());
  std::vector<Eigen::Triplet<typename DerivedV::Scalar>> triplets;
  for (int i = 0; i < (int)vs.size(); ++i)
  {
    triplets.push_back(Eigen::Triplet<typename DerivedV::Scalar>(rs[i],
                                                                 cs[i],
                                                                 vs[i]));
  }
  G.setFromTriplets(triplets.begin(), triplets.end());
  //  std::cout<<"Adjusted"<<fixed<<std::endl;
};

void simplified_covariance_scatter_matrix(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::SparseMatrix<double> &Dx, Eigen::SparseMatrix<double> &Dy,
    Eigen::SparseMatrix<double> &Dz)
{
  using namespace Eigen;
  auto energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
  SparseMatrix<double> Kx, Ky, Kz;
  igl::arap_linear_block(V, F, 0, energy, Kx);
  igl::arap_linear_block(V, F, 1, energy, Ky);
  igl::arap_linear_block(V, F, 2, energy, Kz);
  Dx = Kx.transpose();
  Dy = Ky.transpose();
  Dz = Kz.transpose();
}

void compute_scaffold_gradient_matrix(SCAFData &d_,
                                      Eigen::SparseMatrix<double> &D1,
                                      Eigen::SparseMatrix<double> &D2)
{
  using namespace Eigen;
  Eigen::SparseMatrix<double> G;
  MatrixXi F_s = d_.s_T;
  int vn = d_.v_num;
  MatrixXd V = MatrixXd::Zero(vn, 3);
  V.leftCols(2) = d_.w_uv;
  //  std::cout<<"Avg Mesh Area"<<d_.mesh_measure/d_.mv_num<<std::endl;

  double min_bnd_edge_len = INFINITY;
  int acc_bnd = 0;
  for (int i = 0; i < d_.bnd_sizes.size(); i++)
  {
    int current_size = d_.bnd_sizes[i];

    for (int e = acc_bnd; e < acc_bnd + current_size - 1; e++)
    {
      min_bnd_edge_len = std::min(min_bnd_edge_len,
                                  (d_.w_uv.row(d_.internal_bnd(e)) -
                                   d_.w_uv.row(d_.internal_bnd(e + 1)))
                                      .squaredNorm());
    }
    min_bnd_edge_len = std::min(min_bnd_edge_len,
                                (d_.w_uv.row(d_.internal_bnd(acc_bnd)) -
                                 d_.w_uv.row(d_.internal_bnd(acc_bnd + current_size -
                                                             1)))
                                    .squaredNorm());
    acc_bnd += current_size;
  }

  //  std::cout<<"MinBndEdge"<<min_bnd_edge_len<<std::endl;
  double area_threshold = min_bnd_edge_len / 4.0;

  adjusted_grad(V, F_s, G, area_threshold);
  Eigen::SparseMatrix<double> Dx = G.block(0, 0, F_s.rows(), vn);
  Eigen::SparseMatrix<double> Dy = G.block(F_s.rows(), 0, F_s.rows(), vn);
  Eigen::SparseMatrix<double> Dz = G.block(2 * F_s.rows(), 0, F_s.rows(), vn);

  MatrixXd F1, F2, F3;
  igl::local_basis(V, F_s, F1, F2, F3);
  D1 = F1.col(0).asDiagonal() * Dx + F1.col(1).asDiagonal() * Dy +
       F1.col(2).asDiagonal() * Dz;
  D2 = F2.col(0).asDiagonal() * Dx + F2.col(1).asDiagonal() * Dy +
       F2.col(2).asDiagonal() * Dz;
}


void mesh_improve(igl::SCAFData &s)
{
  using namespace Eigen;
  MatrixXd m_uv = s.w_uv.topRows(s.mv_num);
  MatrixXd V_bnd;
  V_bnd.resize(s.internal_bnd.size(), 2);
  for (int i = 0; i < s.internal_bnd.size(); i++) // redoing step 1.
  {
    V_bnd.row(i) = m_uv.row(s.internal_bnd(i));
  }

  if (s.rect_frame_V.size() == 0)
  {
    Matrix2d ob; // = rect_corners;
    {
      VectorXd uv_max = m_uv.colwise().maxCoeff();
      VectorXd uv_min = m_uv.colwise().minCoeff();
      VectorXd uv_mid = (uv_max + uv_min) / 2.;

      Eigen::Array2d scaf_range(3, 3);
      ob.row(0) = uv_mid.array() + scaf_range * ((uv_min - uv_mid).array());
      ob.row(1) = uv_mid.array() + scaf_range * ((uv_max - uv_mid).array());
    }
    Vector2d rect_len;
    rect_len << ob(1, 0) - ob(0, 0), ob(1, 1) - ob(0, 1);
    int frame_points = 5;

    s.rect_frame_V.resize(4 * frame_points, 2);
    for (int i = 0; i < frame_points; i++)
    {
      // 0,0;0,1
      s.rect_frame_V.row(i) << ob(0, 0), ob(0, 1) + i * rect_len(1) / frame_points;
      // 0,0;1,1
      s.rect_frame_V.row(i + frame_points)
          << ob(0, 0) + i * rect_len(0) / frame_points,
          ob(1, 1);
      // 1,0;1,1
      s.rect_frame_V.row(i + 2 * frame_points) << ob(1, 0), ob(1, 1) - i * rect_len(1) / frame_points;
      // 1,0;0,1
      s.rect_frame_V.row(i + 3 * frame_points)
          << ob(1, 0) - i * rect_len(0) / frame_points,
          ob(0, 1);
      // 0,0;0,1
    }
    s.frame_ids = Eigen::VectorXi::LinSpaced(s.rect_frame_V.rows(), s.mv_num,
                                             s.mv_num +
                                                 s.rect_frame_V.rows());
  }

  // Concatenate Vert and Edge
  MatrixXd V;
  MatrixXi E;
  igl::cat(1, V_bnd, s.rect_frame_V, V);
  E.resize(V.rows(), 2);
  for (int i = 0; i < E.rows(); i++)
    E.row(i) << i, i + 1;
  int acc_bs = 0;
  for (auto bs : s.bnd_sizes)
  {
    E(acc_bs + bs - 1, 1) = acc_bs;
    acc_bs += bs;
  }
  E(V.rows() - 1, 1) = acc_bs;
  assert(acc_bs == s.internal_bnd.size());

  MatrixXd H = MatrixXd::Zero(s.component_sizes.size(), 2);
  {
    int hole_f = 0;
    int hole_i = 0;
    for (auto cs : s.component_sizes)
    {
      for (int i = 0; i < 3; i++)
        H.row(hole_i) += m_uv.row(s.m_T(hole_f, i)); // redoing step 2
      hole_f += cs;
      hole_i++;
    }
  }
  H /= 3.;

  MatrixXd uv2;
  igl::triangle::triangulate(V, E, H, "qYYQ", uv2, s.s_T);
  auto bnd_n = s.internal_bnd.size();

  for (auto i = 0; i < s.s_T.rows(); i++)
    for (auto j = 0; j < s.s_T.cols(); j++)
    {
      auto &x = s.s_T(i, j);
      if (x < bnd_n)
        x = s.internal_bnd(x);
      else
        x += m_uv.rows() - bnd_n;
    }

  igl::cat(1, s.m_T, s.s_T, s.w_T);
  s.w_uv.conservativeResize(m_uv.rows() - bnd_n + uv2.rows(), 2);
  s.w_uv.bottomRows(uv2.rows() - bnd_n) = uv2.bottomRows(-bnd_n + uv2.rows());

  update_scaffold(s);

  auto& d_ = s;
  const auto&v_n = d_.v_num;
  // after_mesh_improve
  if (d_.dim == 2)
  {
    compute_scaffold_gradient_matrix(d_, d_.Dx_s, d_.Dy_s);
  }
  else
  {
    Eigen::SparseMatrix<double> Gs;
    igl::grad(d_.w_uv, d_.s_T, Gs);
    int sf_n = d_.s_T.rows();
    d_.Dx_s = Gs.block(0, 0, sf_n, v_n);
    d_.Dy_s = Gs.block(sf_n, 0, sf_n, v_n);
    d_.Dz_s = Gs.block(2 * sf_n, 0, sf_n, v_n);
  }

  d_.Dx_s.makeCompressed();
  d_.Dy_s.makeCompressed();
  d_.Dz_s.makeCompressed();
  d_.Ri_s = MatrixXd::Zero(d_.Dx_s.rows(),d_.dim * d_.dim);
  d_.Ji_s.resize(d_.Dx_s.rows(), d_.dim * d_.dim);
  d_.W_s.resize(d_.Dx_s.rows(), d_.dim * d_.dim);
}


  void add_new_patch(igl::SCAFData &s, const Eigen::MatrixXd &V_in,
                     const Eigen::MatrixXi &F_ref,
                     const Eigen::RowVectorXd &center)
  {
    using namespace std;
    using namespace Eigen;

    VectorXd M;
    igl::doublearea(V_in, F_ref, M);

    Eigen::MatrixXd V_ref = V_in; // / sqrt(M.sum()/2/igl::PI);
    // M /= M.sum()/igl::PI;
    Eigen::MatrixXd uv_init;
    Eigen::VectorXi bnd;
    Eigen::MatrixXd bnd_uv;

    std::vector<std::vector<int>> all_bnds;
    igl::boundary_loop(F_ref, all_bnds);
    int num_holes = all_bnds.size() - 1;

    std::sort(all_bnds.begin(), all_bnds.end(), [](const std::vector<int> &a, const std::vector<int> &b) {
      return a.size() > b.size();});

    bnd = Map<Eigen::VectorXi>(all_bnds[0].data(), all_bnds[0].size());

    igl::map_vertices_to_circle(V_ref, bnd, bnd_uv);
    bnd_uv *= sqrt(M.sum() / (2 * igl::PI));
    bnd_uv.rowwise() += center;
    s.mesh_measure += M.sum() / 2;
    std::cout << "Mesh Measure" << M.sum() / 2 << std::endl;

    if (num_holes == 0)
    {
      if (bnd.rows() == V_ref.rows())
      {
        std::cout << "All vert on boundary" << std::endl;
        uv_init.resize(V_ref.rows(), 2);
        for (int i = 0; i < bnd.rows(); i++)
        {
          uv_init.row(bnd(i)) = bnd_uv.row(i);
        }
      }
      else
      {
        igl::harmonic(V_ref, F_ref, bnd, bnd_uv, 1, uv_init);
        if (igl::flipped_triangles(uv_init, F_ref).size() != 0)
        {
          std::cout << "Using Uniform Laplacian" << std::endl;
          igl::harmonic(F_ref, bnd, bnd_uv, 1, uv_init); // use uniform laplacian
        }
      }
    }
    else
    {
      auto &F = F_ref;
      auto &V = V_in;
      auto &primary_bnd = bnd;
      // fill holes
      int n_filled_faces = 0;
      int real_F_num = F.rows();
      for (int i = 0; i < num_holes; i++)
      {
        n_filled_faces += all_bnds[i + 1].size();
      }
      MatrixXi F_filled(n_filled_faces + real_F_num, 3);
      F_filled.topRows(real_F_num) = F;

      int new_vert_id = V.rows();
      int new_face_id = real_F_num;

      for (int i = 0; i < num_holes; i++)
      {
        int cur_bnd_size = all_bnds[i + 1].size();
        auto it = all_bnds[i + 1].begin();
        auto back = all_bnds[i + 1].end() - 1;
        F_filled.row(new_face_id++) << *it, *back, new_vert_id;
        while (it != back)
        {
          F_filled.row(new_face_id++)
              << *(it + 1),
              *(it), new_vert_id;
          it++;
        }
        new_vert_id++;
      }
      assert(new_face_id == F_filled.rows());
      assert(new_vert_id == V.rows() + num_holes);

      igl::harmonic(F_filled, primary_bnd, bnd_uv, 1, uv_init);
      uv_init.conservativeResize(V.rows(), 2);
      if (igl::flipped_triangles(uv_init, F_ref).size() != 0)
      {
        std::cout << "Wrong Choice of Outer bnd:" << std::endl;
        //      assert(false&&"Wrong Choice of outer bnd?");
      }
    }

    s.component_sizes.push_back(F_ref.rows());

    MatrixXd m_uv = s.w_uv.topRows(s.mv_num);
    igl::cat(1, m_uv, uv_init, s.w_uv);
    //  s.mv_num =  s.w_uv.rows();

    s.m_M.conservativeResize(s.mf_num + M.size());
    s.m_M.bottomRows(M.size()) = M / 2;

    //  internal_bnd.conservativeResize(internal_bnd.size()+ bnd.size());
    //  internal_bnd.bottomRows(bnd.size()) = bnd.array() + s.mv_num;
    //  bnd_sizes.push_back(bnd.size());

    for (auto cur_bnd : all_bnds)
    {
      s.internal_bnd.conservativeResize(s.internal_bnd.size() + cur_bnd.size());
      s.internal_bnd.bottomRows(cur_bnd.size()) =
          Map<ArrayXi>(cur_bnd.data(), cur_bnd.size()) + s.mv_num;
      s.bnd_sizes.push_back(cur_bnd.size());
    }

    s.m_T.conservativeResize(s.mf_num + F_ref.rows(), 3);
    s.m_T.bottomRows(F_ref.rows()) = F_ref.array() + s.mv_num;
    s.mf_num += F_ref.rows();

    s.m_V.conservativeResize(s.mv_num + V_ref.rows(), 3);
    s.m_V.bottomRows(V_ref.rows()) = V_ref;
    s.mv_num += V_ref.rows();

    s.rect_frame_V = MatrixXd();

    mesh_improve(s);
  }


void add_new_patch(igl::SCAFData &s, const Eigen::MatrixXd &V_in,
  const Eigen::MatrixXi &F_ref,
  const Eigen::RowVectorXd &center,
  const Eigen::MatrixXd& uv_input)
{
  using namespace std;
  using namespace Eigen;

  VectorXd M;
  igl::doublearea(V_in, F_ref, M);

  Eigen::MatrixXd V_ref = V_in; // / sqrt(M.sum()/2/igl::PI);
  // M /= M.sum()/igl::PI;
  Eigen::MatrixXd uv_init;
  Eigen::VectorXi bnd;
  Eigen::MatrixXd bnd_uv;

  std::vector<std::vector<int>> all_bnds;
  igl::boundary_loop(F_ref, all_bnds);
  int num_holes = all_bnds.size() - 1;

  std::sort(all_bnds.begin(), all_bnds.end(), [](const std::vector<int> &a, const std::vector<int> &b) {
  return a.size() > b.size();});

  bnd = Map<Eigen::VectorXi>(all_bnds[0].data(), all_bnds[0].size());

  s.mesh_measure += M.sum() / 2;
  std::cout << "Mesh Measure" << M.sum() / 2 << std::endl;

  uv_init = uv_input;
  if (false)
  {
    if (bnd.rows() == V_ref.rows())
    {
      std::cout << "All vert on boundary" << std::endl;
      uv_init.resize(V_ref.rows(), 2);
      for (int i = 0; i < bnd.rows(); i++)
      {
        uv_init.row(bnd(i)) = bnd_uv.row(i);
      }
    }
    else
    {
      igl::harmonic(V_ref, F_ref, bnd, bnd_uv, 1, uv_init);
      if (igl::flipped_triangles(uv_init, F_ref).size() != 0)
      {
        std::cout << "Using Uniform Laplacian" << std::endl;
        igl::harmonic(F_ref, bnd, bnd_uv, 1, uv_init); // use uniform laplacian
      }
    }
  }
 

  s.component_sizes.push_back(F_ref.rows());

  MatrixXd m_uv = s.w_uv.topRows(s.mv_num);
  igl::cat(1, m_uv, uv_init, s.w_uv);
  //  s.mv_num =  s.w_uv.rows();

  s.m_M.conservativeResize(s.mf_num + M.size());
  s.m_M.bottomRows(M.size()) = M / 2;

  //  internal_bnd.conservativeResize(internal_bnd.size()+ bnd.size());
  //  internal_bnd.bottomRows(bnd.size()) = bnd.array() + s.mv_num;
  //  bnd_sizes.push_back(bnd.size());

  for (auto cur_bnd : all_bnds)
  {
  s.internal_bnd.conservativeResize(s.internal_bnd.size() + cur_bnd.size());
  s.internal_bnd.bottomRows(cur_bnd.size()) =
  Map<ArrayXi>(cur_bnd.data(), cur_bnd.size()) + s.mv_num;
  s.bnd_sizes.push_back(cur_bnd.size());
  }

  s.m_T.conservativeResize(s.mf_num + F_ref.rows(), 3);
  s.m_T.bottomRows(F_ref.rows()) = F_ref.array() + s.mv_num;
  s.mf_num += F_ref.rows();

  s.m_V.conservativeResize(s.mv_num + V_ref.rows(), 3);
  s.m_V.bottomRows(V_ref.rows()) = V_ref;
  s.mv_num += V_ref.rows();
  
  s.rect_frame_V = MatrixXd();

  mesh_improve(s);
 }

void compute_jacobians(SCAFData &d_, const Eigen::MatrixXd &V_new, bool whole)
{
  auto comp_J2 = [](const Eigen::MatrixXd &uv,
                    const Eigen::SparseMatrix<double> &Dx,
                    const Eigen::SparseMatrix<double> &Dy,
                    Eigen::MatrixXd &Ji) {
    // Ji=[D1*u,D2*u,D1*v,D2*v];
    Ji.resize(Dx.rows(), 4);
    Ji.col(0) = Dx * uv.col(0);
    Ji.col(1) = Dy * uv.col(0);
    Ji.col(2) = Dx * uv.col(1);
    Ji.col(3) = Dy * uv.col(1);
  };
  auto comp_J3 = [](const Eigen::MatrixXd &uv,
                    const Eigen::SparseMatrix<double> &Dx,
                    const Eigen::SparseMatrix<double> &Dy,
                    const Eigen::SparseMatrix<double> &Dz,
                    Eigen::MatrixXd &Ji) {
    // Ji=[D1*u,D2*u,D3*u, D1*v,D2*v, D3*v, D1*w,D2*w,D3*w];
    Ji.resize(Dx.rows(), 9);
    Ji.col(0) = Dx * uv.col(0);
    Ji.col(1) = Dy * uv.col(0);
    Ji.col(2) = Dz * uv.col(0);
    Ji.col(3) = Dx * uv.col(1);
    Ji.col(4) = Dy * uv.col(1);
    Ji.col(5) = Dz * uv.col(1);
    Ji.col(6) = Dx * uv.col(2);
    Ji.col(7) = Dy * uv.col(2);
    Ji.col(8) = Dz * uv.col(2);
  };
  Eigen::MatrixXd m_V_new = V_new.topRows(d_.mv_num);
  if (d_.dim == 2)
  {
    comp_J2(m_V_new, d_.Dx_m, d_.Dy_m, d_.Ji_m);
    if (whole)
      comp_J2(V_new, d_.Dx_s, d_.Dy_s, d_.Ji_s);
  }
  else
  {
    comp_J3(m_V_new, d_.Dx_m, d_.Dy_m, d_.Dz_m, d_.Ji_m);
    if (whole)
      comp_J3(V_new, d_.Dx_s, d_.Dy_s, d_.Dz_s, d_.Ji_s);
  }
}

void compute_jacobians(SCAFData &d_, bool whole)
{
  compute_jacobians(d_, d_.w_uv, whole);
}

double compute_energy_from_jacobians(const Eigen::MatrixXd &Ji,
                                     const Eigen::VectorXd &areas,
                                     igl::SLIMData::SLIM_ENERGY energy_type)
{
  double energy = 0;
  int dim = Ji.cols() == 4 ? 2 : 3;
  if (dim == 2)
  {
    Eigen::Matrix<double, 2, 2> ji;
    for (int i = 0; i < Ji.rows(); i++)
    {
      ji(0, 0) = Ji(i, 0);
      ji(0, 1) = Ji(i, 1);
      ji(1, 0) = Ji(i, 2);
      ji(1, 1) = Ji(i, 3);

      typedef Eigen::Matrix<double, 2, 2> Mat2;
      typedef Eigen::Matrix<double, 2, 1> Vec2;
      Mat2 ri, ti, ui, vi;
      Vec2 sing;
      igl::polar_svd(ji, ri, ti, ui, sing, vi);
      double s1 = sing(0);
      double s2 = sing(1);

      switch (energy_type)
      {
      case igl::SLIMData::ARAP:
      {
        energy += areas(i) * (pow(s1 - 1, 2) + pow(s2 - 1, 2));
        break;
      }
      case igl::SLIMData::SYMMETRIC_DIRICHLET:
      {
        energy +=
            areas(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) - 4);
        break;
      }
      case igl::SLIMData::LOG_ARAP:
      {
        energy += areas(i) * (pow(log(s1), 2) + pow(log(s2), 2));
        break;
      }
      case igl::SLIMData::CONFORMAL:
      {
        energy += areas(i) * ((pow(s1, 2) + pow(s2, 2)) / (2 * s1 * s2));
        break;
      }
      default:
        break;
      }
    }
  }
  else
  {
    Eigen::Matrix<double, 3, 3> ji;
    for (int i = 0; i < Ji.rows(); i++)
    {
      ji(0, 0) = Ji(i, 0);
      ji(0, 1) = Ji(i, 1);
      ji(0, 2) = Ji(i, 2);
      ji(1, 0) = Ji(i, 3);
      ji(1, 1) = Ji(i, 4);
      ji(1, 2) = Ji(i, 5);
      ji(2, 0) = Ji(i, 6);
      ji(2, 1) = Ji(i, 7);
      ji(2, 2) = Ji(i, 8);

      typedef Eigen::Matrix<double, 3, 3> Mat3;
      typedef Eigen::Matrix<double, 3, 1> Vec3;
      Mat3 ri, ti, ui, vi;
      Vec3 sing;
      igl::polar_svd(ji, ri, ti, ui, sing, vi);
      double s1 = sing(0);
      double s2 = sing(1);
      double s3 = sing(2);

      switch (energy_type)
      {
      case igl::SLIMData::ARAP:
      {
        energy += areas(i) * (pow(s1 - 1, 2) + pow(s2 - 1, 2) + pow(s3 - 1, 2));
        break;
      }
      case igl::SLIMData::SYMMETRIC_DIRICHLET:
      {
        energy += areas(i) * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2) + pow(s3, 2) + pow(s3, -2));
        break;
      }
      case igl::SLIMData::LOG_ARAP:
      {
        energy += areas(i) * (pow(log(s1), 2) + pow(log(std::abs(s2)), 2) + pow(log(std::abs(s3)), 2));
        break;
      }
      case igl::SLIMData::CONFORMAL:
      {
        energy += areas(i) * ((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow(s1 * s2 * s3, 2. / 3.)));
        break;
      }
      case igl::SLIMData::EXP_CONFORMAL:
      {
        energy += areas(i) * exp((pow(s1, 2) + pow(s2, 2) + pow(s3, 2)) / (3 * pow(s1 * s2 * s3, 2. / 3.)));
        break;
      }
      default:
        assert(false);
      }
    }
  }
  return energy;
}

double compute_soft_constraint_energy(const SCAFData &d_)
{
  double e = 0;
  for (auto const &x : d_.soft_cons)
    e += d_.soft_const_p * (x.second - d_.w_uv.row(x.first)).squaredNorm();

  return e;
}

double compute_energy(SCAFData &d_, bool whole)
{
  compute_jacobians(d_, whole);
  double energy = compute_energy_from_jacobians(d_.Ji_m, d_.m_M, d_.slim_energy);

  if (whole)
    energy += compute_energy_from_jacobians(d_.Ji_s, d_.s_M, d_.scaf_energy);
  energy += compute_soft_constraint_energy(d_);
  return energy;
}

double compute_energy(SCAFData &d_, Eigen::MatrixXd &w_uv)
{
  bool whole = w_uv.rows() == d_.v_num;
  compute_jacobians(d_, w_uv, whole);
  double energy = compute_energy_from_jacobians(d_.Ji_m, d_.m_M, d_.slim_energy);

  if (whole)
    energy += compute_energy_from_jacobians(d_.Ji_s, d_.s_M, d_.scaf_energy);
  energy += compute_soft_constraint_energy(d_);
  return energy;
}

void adjust_scaf_weight(SCAFData &d_, double new_weight)
{
  d_.scaffold_factor = new_weight;
  update_scaffold(d_);
}

template <int dim>
void update_weights_and_closest_rotations(
    const Eigen::MatrixXd &Ji,
    igl::SLIMData::SLIM_ENERGY energy_type,
    Eigen::MatrixXd &W,
    Eigen::MatrixXd &Ri);

template <>
void update_weights_and_closest_rotations<2>(
    const Eigen::MatrixXd &Ji,
    igl::SLIMData::SLIM_ENERGY energy_type,
    Eigen::MatrixXd &W,
    Eigen::MatrixXd &Ri)
{
  const double eps = 1e-8;

  W.resize(Ji.rows(), 4);
  Ri.resize(Ji.rows(), 4);
  for (int i = 0; i < Ji.rows(); ++i)
  {
    typedef Eigen::Matrix<double, 2, 2> Mat2;
    typedef Eigen::Matrix<double, 2, 1> Vec2;
    Mat2 ji, ri, ti, ui, vi;
    Vec2 sing;
    Vec2 closest_sing_vec;
    Mat2 mat_W;
    Vec2 m_sing_new;
    double s1, s2;

    ji(0, 0) = Ji(i, 0);
    ji(0, 1) = Ji(i, 1);
    ji(1, 0) = Ji(i, 2);
    ji(1, 1) = Ji(i, 3);

    igl::polar_svd(ji, ri, ti, ui, sing, vi);

    s1 = sing(0);
    s2 = sing(1);

    // Branch between mesh face and scaffold faces.
    switch (energy_type)
    {
    case igl::SLIMData::ARAP:
    {
      m_sing_new << 1, 1;
      break;
    }
    case igl::SLIMData::SYMMETRIC_DIRICHLET:
    {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(
                                                     s2_g / (2 * (s2 - 1)));
      break;
    }
    case igl::SLIMData::LOG_ARAP:
    {
      double s1_g = 2 * (log(s1) / s1);
      double s2_g = 2 * (log(s2) / s2);
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(
                                                     s2_g / (2 * (s2 - 1)));
      break;
    }
    case igl::SLIMData::CONFORMAL:
    {
      double s1_g = 1 / (2 * s2) - s2 / (2 * pow(s1, 2));
      double s2_g = 1 / (2 * s1) - s1 / (2 * pow(s2, 2));

      double geo_avg = sqrt(s1 * s2);
      double s1_min = geo_avg;
      double s2_min = geo_avg;

      m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(
                                                          s2_g / (2 * (s2 - s2_min)));

      // change local step
      closest_sing_vec << s1_min, s2_min;
      ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
      break;
    }
    default:
      break;
    }

    if (abs(s1 - 1) < eps)
      m_sing_new(0) = 1;
    if (abs(s2 - 1) < eps)
      m_sing_new(1) = 1;

    mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

    W(i, 0) = mat_W(0, 0);
    W(i, 1) = mat_W(0, 1);
    W(i, 2) = mat_W(1, 0);
    W(i, 3) = mat_W(1, 1);

    // 2) Update local step (doesn't have to be a rotation, for instance in case of conformal energy)
    Ri(i, 0) = ri(0, 0);
    Ri(i, 1) = ri(1, 0);
    Ri(i, 2) = ri(0, 1);
    Ri(i, 3) = ri(1, 1);
  }
}

template <>
void update_weights_and_closest_rotations<3>(
    const Eigen::MatrixXd &Ji,
    igl::SLIMData::SLIM_ENERGY energy_type,
    Eigen::MatrixXd &W,
    Eigen::MatrixXd &Ri)
{
  const double eps = 1e-8;

  typedef Eigen::Matrix<double, 3, 1> Vec3;
  typedef Eigen::Matrix<double, 3, 3> Mat3;
  Mat3 ji;
  Vec3 m_sing_new;
  Vec3 closest_sing_vec;
  const double sqrt_2 = sqrt(2);

  W.resize(Ji.rows(), 9);
  Ri.resize(Ji.rows(), 9);
  for (int i = 0; i < Ji.rows(); ++i)
  {
    ji(0, 0) = Ji(i, 0);
    ji(0, 1) = Ji(i, 1);
    ji(0, 2) = Ji(i, 2);
    ji(1, 0) = Ji(i, 3);
    ji(1, 1) = Ji(i, 4);
    ji(1, 2) = Ji(i, 5);
    ji(2, 0) = Ji(i, 6);
    ji(2, 1) = Ji(i, 7);
    ji(2, 2) = Ji(i, 8);

    Mat3 ri, ti, ui, vi;
    Vec3 sing;
    igl::polar_svd(ji, ri, ti, ui, sing, vi);

    double s1 = sing(0);
    double s2 = sing(1);
    double s3 = sing(2);

    // 1) Update Weights
    switch (energy_type)
    {
    case igl::SLIMData::ARAP:
    {
      m_sing_new << 1, 1, 1;
      break;
    }
    case igl::SLIMData::LOG_ARAP:
    {
      double s1_g = 2 * (log(s1) / s1);
      double s2_g = 2 * (log(s2) / s2);
      double s3_g = 2 * (log(s3) / s3);
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(
                                                     s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));
      break;
    }
    case igl::SLIMData::SYMMETRIC_DIRICHLET:
    {
      double s1_g = 2 * (s1 - pow(s1, -3));
      double s2_g = 2 * (s2 - pow(s2, -3));
      double s3_g = 2 * (s3 - pow(s3, -3));
      m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(
                                                     s2_g / (2 * (s2 - 1))),
          sqrt(s3_g / (2 * (s3 - 1)));
      break;
    }

    case igl::SLIMData::CONFORMAL:
    {
      double common_div = 9 * (pow(s1 * s2 * s3, 5. / 3.));

      double s1_g =
          (-2 * s2 * s3 * (pow(s2, 2) + pow(s3, 2) - 2 * pow(s1, 2))) / common_div;
      double s2_g =
          (-2 * s1 * s3 * (pow(s1, 2) + pow(s3, 2) - 2 * pow(s2, 2))) / common_div;
      double s3_g =
          (-2 * s1 * s2 * (pow(s1, 2) + pow(s2, 2) - 2 * pow(s3, 2))) / common_div;

      double closest_s = sqrt(pow(s1, 2) + pow(s3, 2)) / sqrt_2;
      double s1_min = closest_s;
      double s2_min = closest_s;
      double s3_min = closest_s;

      m_sing_new << sqrt(s1_g / (2 * (s1 - s1_min))), sqrt(
                                                          s2_g / (2 * (s2 - s2_min))),
          sqrt(
              s3_g / (2 * (s3 - s3_min)));

      // change local step
      closest_sing_vec << s1_min, s2_min, s3_min;
      ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
      break;
    }

    default:
      assert(false);
    }
    if (std::abs(s1 - 1) < eps)
      m_sing_new(0) = 1;
    if (std::abs(s2 - 1) < eps)
      m_sing_new(1) = 1;
    if (std::abs(s3 - 1) < eps)
      m_sing_new(2) = 1;
    Mat3 mat_W;
    mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();

    W(i, 0) = mat_W(0, 0);
    W(i, 1) = mat_W(0, 1);
    W(i, 2) = mat_W(0, 2);
    W(i, 3) = mat_W(1, 0);
    W(i, 4) = mat_W(1, 1);
    W(i, 5) = mat_W(1, 2);
    W(i, 6) = mat_W(2, 0);
    W(i, 7) = mat_W(2, 1);
    W(i, 8) = mat_W(2, 2);

    // 2) Update closest rotations (not rotations in case of conformal energy)
    Ri(i, 0) = ri(0, 0);
    Ri(i, 1) = ri(1, 0);
    Ri(i, 2) = ri(2, 0);
    Ri(i, 3) = ri(0, 1);
    Ri(i, 4) = ri(1, 1);
    Ri(i, 5) = ri(2, 1);
    Ri(i, 6) = ri(0, 2);
    Ri(i, 7) = ri(1, 2);
    Ri(i, 8) = ri(2, 2);
  }
}

void buildAm(const Eigen::VectorXd &sqrt_M,
             const Eigen::SparseMatrix<double> &Dx,
             const Eigen::SparseMatrix<double> &Dy,
             const Eigen::MatrixXd &W,
             Eigen::SparseMatrix<double> &Am)
{
  using namespace Eigen;
  // formula (35) in paper
  std::vector<Triplet<double>> IJV;
  const int f_n = W.rows();
  const int v_n = Dx.cols();

  IJV.reserve(4 * (Dx.outerSize() + Dy.outerSize()));

  /*A = [W11*Dx, W12*Dx;
  W11*Dy, W12*Dy;
  W21*Dx, W22*Dx;
  W21*Dy, W22*Dy];*/
  for (int k = 0; k < Dx.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it)
    {
      int dx_r = it.row();
      int dx_c = it.col();
      double val = it.value() * sqrt_M(dx_r);

      IJV.push_back(Triplet<double>(dx_r, dx_c, val * W(dx_r, 0)));
      IJV.push_back(Triplet<double>(dx_r, v_n + dx_c, val * W(dx_r, 1)));

      IJV.push_back(Triplet<double>(2 * f_n + dx_r, dx_c, val * W(dx_r, 2)));
      IJV.push_back(
          Triplet<double>(2 * f_n + dx_r, v_n + dx_c, val * W(dx_r, 3)));
    }
  }

  for (int k = 0; k < Dy.outerSize(); ++k)
  {
    for (SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it)
    {
      int dy_r = it.row();
      int dy_c = it.col();
      double val = it.value() * sqrt_M(dy_r);

      IJV.push_back(Triplet<double>(f_n + dy_r, dy_c,
                                    val * W(dy_r, 0)));
      IJV.push_back(Triplet<double>(f_n + dy_r, v_n + dy_c,
                                    val * W(dy_r, 1)));

      IJV.push_back(Triplet<double>(3 * f_n + dy_r, dy_c,
                                    val * W(dy_r, 2)));
      IJV.push_back(Triplet<double>(3 * f_n + dy_r, v_n + dy_c,
                                    val * W(dy_r, 3)));
    }
  }
  Am.setFromTriplets(IJV.begin(), IJV.end());
}

void buildAm(const Eigen::VectorXd &sqrt_M,
             const Eigen::SparseMatrix<double> &Dx,
             const Eigen::SparseMatrix<double> &Dy,
             const Eigen::SparseMatrix<double> &Dz,
             const Eigen::MatrixXd &W,
             Eigen::SparseMatrix<double> &Am)
{
  using namespace Eigen;
  // formula (35) in paper
  std::vector<Triplet<double>> IJV;
  IJV.reserve(9 * (Dx.outerSize() + Dy.outerSize() + Dz.outerSize()));

  const int f_n = W.rows();
  const int v_n = Dx.cols();

  /*A = [W11*Dx, W12*Dx, W13*Dx;
         W11*Dy, W12*Dy, W13*Dy;
         W11*Dz, W12*Dz, W13*Dz;
         W21*Dx, W22*Dx, W23*Dx;
         W21*Dy, W22*Dy, W23*Dy;
         W21*Dz, W22*Dz, W23*Dz;
         W31*Dx, W32*Dx, W33*Dx;
         W31*Dy, W32*Dy, W33*Dy;
         W31*Dz, W32*Dz, W33*Dz;];*/
  for (int k = 0; k < Dx.outerSize(); k++)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dx, k); it; ++it)
    {
      int dx_r = it.row();
      int dx_c = it.col();
      double val = it.value();

      double m_0 = sqrt_M(dx_r);
      double m_3 = sqrt_M(dx_r);
      double m_6 = sqrt_M(dx_r);
      IJV.push_back(Eigen::Triplet<double>(dx_r, dx_c, m_0 * val * W(dx_r, 0)));
      IJV.push_back(Eigen::Triplet<double>(dx_r,
                                           v_n + dx_c,
                                           m_0 * val * W(dx_r, 1)));
      IJV.push_back(Eigen::Triplet<double>(dx_r,
                                           2 * v_n + dx_c,
                                           m_0 * val * W(dx_r, 2)));

      IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r,
                                           dx_c,
                                           m_3 * val * W(dx_r, 3)));
      IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r,
                                           v_n + dx_c,
                                           m_3 * val * W(dx_r, 4)));
      IJV.push_back(Eigen::Triplet<double>(3 * f_n + dx_r,
                                           2 * v_n + dx_c,
                                           m_3 * val * W(dx_r, 5)));

      IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r,
                                           dx_c,
                                           m_6 * val * W(dx_r, 6)));
      IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r,
                                           v_n + dx_c,
                                           m_6 * val * W(dx_r, 7)));
      IJV.push_back(Eigen::Triplet<double>(6 * f_n + dx_r,
                                           2 * v_n + dx_c,
                                           m_6 * val * W(dx_r, 8)));
    }
  }

  for (int k = 0; k < Dy.outerSize(); k++)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dy, k); it; ++it)
    {
      int dy_r = it.row();
      int dy_c = it.col();
      double val = it.value();

      double m_1 = sqrt_M(dy_r);
      double m_4 = sqrt_M(dy_r);
      double m_7 = sqrt_M(dy_r);
      IJV.push_back(Eigen::Triplet<double>(f_n + dy_r,
                                           dy_c,
                                           m_1 * val * W(dy_r, 0)));
      IJV.push_back(Eigen::Triplet<double>(f_n + dy_r,
                                           v_n + dy_c,
                                           m_1 * val * W(dy_r, 1)));
      IJV.push_back(Eigen::Triplet<double>(f_n + dy_r,
                                           2 * v_n + dy_c,
                                           m_1 * val * W(dy_r, 2)));

      IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r,
                                           dy_c,
                                           m_4 * val * W(dy_r, 3)));
      IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r,
                                           v_n + dy_c,
                                           m_4 * val * W(dy_r, 4)));
      IJV.push_back(Eigen::Triplet<double>(4 * f_n + dy_r,
                                           2 * v_n + dy_c,
                                           m_4 * val * W(dy_r, 5)));

      IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r,
                                           dy_c,
                                           m_7 * val * W(dy_r, 6)));
      IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r,
                                           v_n + dy_c,
                                           m_7 * val * W(dy_r, 7)));
      IJV.push_back(Eigen::Triplet<double>(7 * f_n + dy_r,
                                           2 * v_n + dy_c,
                                           m_7 * val * W(dy_r, 8)));
    }
  }

  for (int k = 0; k < Dz.outerSize(); k++)
  {
    for (Eigen::SparseMatrix<double>::InnerIterator it(Dz, k); it; ++it)
    {
      int dz_r = it.row();
      int dz_c = it.col();
      double val = it.value();
      double m_2 = sqrt_M(dz_r);
      double m_5 = sqrt_M(dz_r);
      double m_8 = sqrt_M(dz_r);
      IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r,
                                           dz_c,
                                           m_2 * val * W(dz_r, 0)));
      IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r,
                                           v_n + dz_c,
                                           m_2 * val * W(dz_r, 1)));
      IJV.push_back(Eigen::Triplet<double>(2 * f_n + dz_r,
                                           2 * v_n + dz_c,
                                           m_2 * val * W(dz_r, 2)));

      IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r,
                                           dz_c,
                                           m_5 * val * W(dz_r, 3)));
      IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r,
                                           v_n + dz_c,
                                           m_5 * val * W(dz_r, 4)));
      IJV.push_back(Eigen::Triplet<double>(5 * f_n + dz_r,
                                           2 * v_n + dz_c,
                                           m_5 * val * W(dz_r, 5)));

      IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r,
                                           dz_c,
                                           m_8 * val * W(dz_r, 6)));
      IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r,
                                           v_n + dz_c,
                                           m_8 * val * W(dz_r, 7)));
      IJV.push_back(Eigen::Triplet<double>(8 * f_n + dz_r,
                                           2 * v_n + dz_c,
                                           m_8 * val * W(dz_r, 8)));
    }
  }

  Am.setFromTriplets(IJV.begin(), IJV.end());
}

void buildRhs(const Eigen::VectorXd &sqrt_M,
              const Eigen::MatrixXd &W,
              const Eigen::MatrixXd &Ri,
              Eigen::VectorXd &f_rhs)
{
  const int dim = (W.cols() == 4) ? 2 : 3;
  const int f_n = W.rows();
  f_rhs.resize(dim * dim * f_n);

  if (dim == 2)
  {
    /*b = [W11*R11 + W12*R21; (formula (36))
           W11*R12 + W12*R22;
           W21*R11 + W22*R21;
           W21*R12 + W22*R22];*/
    for (int i = 0; i < f_n; i++)
    {
      auto sqrt_area = sqrt_M(i);
      f_rhs(i + 0 * f_n) = sqrt_area * (W(i, 0) * Ri(i, 0) + W(i, 1) * Ri(i, 1));
      f_rhs(i + 1 * f_n) = sqrt_area * (W(i, 0) * Ri(i, 2) + W(i, 1) * Ri(i, 3));
      f_rhs(i + 2 * f_n) = sqrt_area * (W(i, 2) * Ri(i, 0) + W(i, 3) * Ri(i, 1));
      f_rhs(i + 3 * f_n) = sqrt_area * (W(i, 2) * Ri(i, 2) + W(i, 3) * Ri(i, 3));
    }
  }
  else
  {
    /*b = [W11*R11 + W12*R21 + W13*R31;
           W11*R12 + W12*R22 + W13*R32;
           W11*R13 + W12*R23 + W13*R33;
           W21*R11 + W22*R21 + W23*R31;
           W21*R12 + W22*R22 + W23*R32;
           W21*R13 + W22*R23 + W23*R33;
           W31*R11 + W32*R21 + W33*R31;
           W31*R12 + W32*R22 + W33*R32;
           W31*R13 + W32*R23 + W33*R33;];*/
    for (int i = 0; i < f_n; i++)
    {
      auto sqrt_area = sqrt_M(i);
      f_rhs(i + 0 * f_n) = sqrt_area *
                           (W(i, 0) * Ri(i, 0) + W(i, 1) * Ri(i, 1) + W(i, 2) * Ri(i, 2));
      f_rhs(i + 1 * f_n) = sqrt_area *
                           (W(i, 0) * Ri(i, 3) + W(i, 1) * Ri(i, 4) + W(i, 2) * Ri(i, 5));
      f_rhs(i + 2 * f_n) = sqrt_area *
                           (W(i, 0) * Ri(i, 6) + W(i, 1) * Ri(i, 7) + W(i, 2) * Ri(i, 8));
      f_rhs(i + 3 * f_n) = sqrt_area *
                           (W(i, 3) * Ri(i, 0) + W(i, 4) * Ri(i, 1) + W(i, 5) * Ri(i, 2));
      f_rhs(i + 4 * f_n) = sqrt_area *
                           (W(i, 3) * Ri(i, 3) + W(i, 4) * Ri(i, 4) + W(i, 5) * Ri(i, 5));
      f_rhs(i + 5 * f_n) = sqrt_area *
                           (W(i, 3) * Ri(i, 6) + W(i, 4) * Ri(i, 7) + W(i, 5) * Ri(i, 8));
      f_rhs(i + 6 * f_n) = sqrt_area *
                           (W(i, 6) * Ri(i, 0) + W(i, 7) * Ri(i, 1) + W(i, 8) * Ri(i, 2));
      f_rhs(i + 7 * f_n) = sqrt_area *
                           (W(i, 6) * Ri(i, 3) + W(i, 7) * Ri(i, 4) + W(i, 8) * Ri(i, 5));
      f_rhs(i + 8 * f_n) = sqrt_area *
                           (W(i, 6) * Ri(i, 6) + W(i, 7) * Ri(i, 7) + W(i, 8) * Ri(i, 8));
    }
  }
}

void sparse_slice(Eigen::SparseMatrix<double> & A,
                  const Eigen::VectorXi & unknown_ids,
                  const Eigen::VectorXi &known_ids,
                  Eigen::SparseMatrix<double>& Au, Eigen::SparseMatrix<double> & Ae)
{
  using namespace Eigen;
  using TY = double;
  using TX = double;
  auto &X = A;

  int xm = X.rows();
  int xn = X.cols();
  int ym = xm;
  int yn = unknown_ids.size();
  int ykn = known_ids.size();

  std::vector<int> CI(xn, -1);
  std::vector<int> CKI(xn, -1);
  // initialize to -1
  for (int i = 0; i < yn; i++)
    CI[unknown_ids(i)] = (i);
  for (int i = 0; i < ykn; i++)
    CKI[known_ids(i)] = i;
  Eigen::DynamicSparseMatrix<TY, Eigen::ColMajor> dyn_Y(ym, yn);
  Eigen::DynamicSparseMatrix<TY, Eigen::ColMajor> dyn_K(ym, ykn);
  // Take a guess at the number of nonzeros (this assumes uniform distribution
  // not banded or heavily diagonal)
  dyn_Y.reserve(A.nonZeros());
  dyn_K.reserve(A.nonZeros() * ykn / xn);
  // Iterate over outside
  for (int k = 0; k < X.outerSize(); ++k)
  {
    // Iterate over inside
    if (CI[k] != -1)
      for (typename Eigen::SparseMatrix<double>::InnerIterator it(X, k); it;
           ++it)
      {
        dyn_Y.coeffRef(it.row(), CI[it.col()]) = it.value();
      }
    else
      for (typename Eigen::SparseMatrix<TX>::InnerIterator it(X, k); it;
           ++it)
      {
        dyn_K.coeffRef(it.row(), CKI[it.col()]) = it.value();
      }
  }
  Au = Eigen::SparseMatrix<TY>(dyn_Y);
  Ae = Eigen::SparseMatrix<double>(dyn_K);
}

void get_complement(const Eigen::VectorXi& bnd_ids, int v_n, Eigen::ArrayXi& unknown_ids)
{ // get the complement of bnd_ids.
  int assign = 0, i = 0;
  for (int get = 0; i < v_n && get < bnd_ids.size(); i++)
{
  if (bnd_ids(get) == i)
  get++;
  else
  unknown_ids(assign++) = i;
}
  while (i < v_n)
  unknown_ids(assign++) = i++;
  assert(assign + bnd_ids.size() == v_n);
}


void build_surface_linear_system(const SCAFData &d_, Eigen::SparseMatrix<double> &L, Eigen::VectorXd &rhs)
{

  using namespace Eigen;
  using namespace std;

  const int v_n = d_.v_num - (d_.frame_ids.size());
  const int dim = d_.dim;
  const int f_n = d_.mf_num;

  // to get the  complete A
  Eigen::VectorXd sqrtM = d_.m_M.array().sqrt();
  Eigen::SparseMatrix<double> A(dim * dim * f_n, dim * v_n);
  auto decoy_Dx_m = d_.Dx_m;
  decoy_Dx_m.conservativeResize(d_.W_m.rows(), v_n);
  auto decoy_Dy_m = d_.Dy_m;
  decoy_Dy_m.conservativeResize(d_.W_m.rows(), v_n);
  if (dim == 2)
  {
    buildAm(sqrtM, decoy_Dx_m, decoy_Dy_m, d_.W_m, A);
  }
  else
  {
    auto decoy_Dz_m = d_.Dz_m;
    decoy_Dz_m.conservativeResize(d_.W_m.rows(), v_n);
    buildAm(sqrtM, decoy_Dx_m, decoy_Dy_m, decoy_Dz_m, d_.W_m, A);
  }

  const VectorXi & bnd_ids = d_.fixed_ids;
  auto bnd_n = bnd_ids.size();
  if (bnd_n == 0) {

    Eigen::SparseMatrix<double> At = A.transpose();
    At.makeCompressed();

    Eigen::SparseMatrix<double> id_m(At.rows(), At.rows());
    id_m.setIdentity();

    L = At * A;

    Eigen::VectorXd frhs;
    buildRhs(sqrtM, d_.W_m, d_.Ri_m, frhs);
    rhs = At * frhs;
  } else {
    MatrixXd bnd_pos;
    igl::slice(d_.w_uv, bnd_ids, 1, bnd_pos);
    ArrayXi known_ids(bnd_ids.size() * dim);
    ArrayXi unknown_ids((v_n - bnd_ids.rows()) * dim);
    get_complement(bnd_ids, v_n, unknown_ids);
    VectorXd known_pos(bnd_ids.size() * dim);
    for (int d = 0; d < dim; d++)
    {
      auto n_b = bnd_ids.rows();
      known_ids.segment(d * n_b, n_b) = bnd_ids.array() + d * v_n;
      known_pos.segment(d * n_b, n_b) = bnd_pos.col(d);
      unknown_ids.block(d * (v_n - n_b), 0, v_n - n_b, unknown_ids.cols()) =
          unknown_ids.topRows(v_n - n_b) + d * v_n;
    }

    Eigen::SparseMatrix<double> Au, Ae;
    sparse_slice(A, unknown_ids, known_ids, Au, Ae);

    Eigen::SparseMatrix<double> Aut = Au.transpose();
    Aut.makeCompressed();

    L = Aut * Au;

    Eigen::VectorXd frhs;
    buildRhs(sqrtM, d_.W_m, d_.Ri_m, frhs);

    rhs = Aut * (frhs - Ae * known_pos);
  }

  // add soft constraints.
  for (auto const &x : d_.soft_cons)
  {
    int v_idx = x.first;

    for (int d = 0; d < dim; d++)
    {
      rhs(d * (v_n) + v_idx) += d_.soft_const_p * x.second(d); // rhs
      L.coeffRef(d * v_n + v_idx,
                 d * v_n + v_idx) += d_.soft_const_p; // diagonal
    }
  }
}

void build_scaffold_linear_system(const SCAFData &d_, Eigen::SparseMatrix<double> &L, Eigen::VectorXd &rhs)
{
  using namespace Eigen;

  const int f_n = d_.W_s.rows();
  const int v_n = d_.Dx_s.cols();
  const int dim = d_.dim;

  Eigen::VectorXd sqrtM = d_.s_M.array().sqrt();
  Eigen::SparseMatrix<double> A(dim * dim * f_n, dim * v_n);
  if (dim == 2)
    buildAm(sqrtM, d_.Dx_s, d_.Dy_s, d_.W_s, A);
  else
    buildAm(sqrtM, d_.Dx_s, d_.Dy_s, d_.Dz_s, d_.W_s, A);

//  const VectorXi &bnd_ids = d_.frame_ids;
  VectorXi bnd_ids;
  igl::cat(1, d_.fixed_ids, d_.frame_ids, bnd_ids);

  auto bnd_n = bnd_ids.size();
  assert(bnd_n > 0);
  MatrixXd bnd_pos;
  igl::slice(d_.w_uv, bnd_ids, 1, bnd_pos);

  ArrayXi known_ids(bnd_ids.size() * dim);
  ArrayXi unknown_ids((v_n - bnd_ids.rows()) * dim);

  get_complement(bnd_ids, v_n, unknown_ids);

  VectorXd known_pos(bnd_ids.size() * dim);
  for (int d = 0; d < dim; d++)
  {
    auto n_b = bnd_ids.rows();
    known_ids.segment(d * n_b, n_b) = bnd_ids.array() + d * v_n;
    known_pos.segment(d * n_b, n_b) = bnd_pos.col(d);
    unknown_ids.block(d * (v_n - n_b), 0, v_n - n_b, unknown_ids.cols()) =
        unknown_ids.topRows(v_n - n_b) + d * v_n;
  }
  Eigen::VectorXd sqrt_M = d_.s_M.array().sqrt();

  // slice
  // 'manual slicing for A(:, unknown/known)'
  Eigen::SparseMatrix<double> Au, Ae;
  sparse_slice(A, unknown_ids, known_ids, Au, Ae);

  Eigen::SparseMatrix<double> Aut = Au.transpose();
  Aut.makeCompressed();

  L = Aut * Au;

  Eigen::VectorXd frhs;
  buildRhs(sqrtM, d_.W_s, d_.Ri_s, frhs);

  rhs = Aut * (frhs - Ae * known_pos);
}

void solve_weighted_arap(SCAFData &d_, Eigen::MatrixXd &uv)
{
  using namespace Eigen;
  using namespace std;
  int dim = d_.dim;
  igl::Timer timer;
  timer.start();

  VectorXi bnd_ids;
  igl::cat(1, d_.fixed_ids, d_.frame_ids, bnd_ids);
  const auto v_n = d_.v_num;
  const auto bnd_n = bnd_ids.size();
  assert(bnd_n > 0);
  MatrixXd bnd_pos;
  igl::slice(d_.w_uv, bnd_ids, 1, bnd_pos);

  ArrayXi known_ids(bnd_n * dim);
  ArrayXi unknown_ids((v_n - bnd_n) * dim);
//  std::cout<<bnd_ids<<std::endl;

  get_complement(bnd_ids, v_n, unknown_ids);

  VectorXd known_pos(bnd_ids.size() * dim);
  for (int d = 0; d < dim; d++)
  {
    auto n_b = bnd_ids.rows();
    known_ids.segment(d * n_b, n_b) = bnd_ids.array() + d * v_n;
    known_pos.segment(d * n_b, n_b) = bnd_pos.col(d);
    unknown_ids.block(d * (v_n - n_b), 0, v_n - n_b, unknown_ids.cols()) =
        unknown_ids.topRows(v_n - n_b) + d * v_n;
  }
  //std::cout<<"Slicing Knowns "<<timer.getElapsedTime()<<std::endl;
  //timer.start();

  Eigen::SparseMatrix<double> L;
  Eigen::VectorXd rhs;

  // fixed frame solving:
  // x_e as the fixed frame, x_u for unknowns (mesh + unknown scaffold)
  // min ||(A_u*x_u + A_e*x_e) - b||^2
  // => A_u'*A_u*x_u + A_u'*A_e*x_e = Au'*b
  // => A_u'*A_u*x_u  = Au'* (b - A_e*x_e) := Au'* b_u
  // => L * x_u = rhs
  //
  // separate matrix build:
  // min ||A_m x_m - b_m||^2 + ||A_s x_all - b_s||^2 + soft + proximal
  // First change dimension of A_m to fit for x_all
  // (Not just at the end, since x_all is flattened along dimensions)
  // L = A_m'*A_m + A_s'*A_s + soft + proximal
  // rhs = A_m'* b_m + A_s' * b_s + soft + proximal
  //
  using namespace std;
  Eigen::SparseMatrix<double> L_m, L_s;
  Eigen::VectorXd rhs_m, rhs_s;
  build_surface_linear_system(d_, L_m, rhs_m);  // complete Am, with soft
  build_scaffold_linear_system(d_, L_s, rhs_s); // complete As, without proximal
  // we don't need proximal term

  L = L_m + L_s;
  rhs = rhs_m + rhs_s;
  L.makeCompressed();

  Eigen::VectorXd unknown_Uc((v_n - d_.frame_ids.size() - d_.fixed_ids.size()) * dim), Uc(dim * v_n);
  bool solve_with_cg = (d_.dim == 3); // use CG in 3D
  if (solve_with_cg)
  {
    for (auto t : {1e-6})
    {
      timer.start();
      ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Upper>
          CGsolver;
      CGsolver.setTolerance(t);
      unknown_Uc = CGsolver.compute(L).solve(rhs);
      cout << t << "CGSolve = " << timer.getElapsedTime() << endl;
      std::cout << "#iterations:     " << CGsolver.iterations() << std::endl;
      std::cout << "estimated error: " << CGsolver.error() << std::endl;
    }
  }
  else
  {
    SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    unknown_Uc = solver.compute(L).solve(rhs);
    //cout << "Direct Solve = " << timer.getElapsedTime() << endl;
  }
  //timer.start();
  igl::slice_into(unknown_Uc, unknown_ids.matrix(), 1, Uc);
  igl::slice_into(known_pos, known_ids.matrix(), 1, Uc);

  uv = Map<Matrix<double, -1, -1, Eigen::ColMajor>>(Uc.data(), v_n, dim);
  //for (int i = 0; i < dim; i++)
  //  uv.col(i) = Uc.block(i * v_n, 0, v_n, 1);
  //cout << "Slice back = " << timer.getElapsedTime() << endl;
}

double perform_iteration(SCAFData &d_)
{
  auto &w_uv = d_.w_uv;
  Eigen::MatrixXd V_out = w_uv;
  //solve_weighted_proxy(V_out);
  {
    auto &uv_new = V_out;
    compute_jacobians(d_, uv_new, true);
    if (d_.dim == 2)
    {
      update_weights_and_closest_rotations<2>(d_.Ji_m, d_.slim_energy,
                                              d_.W_m, d_.Ri_m);
      update_weights_and_closest_rotations<2>(d_.Ji_s, d_.scaf_energy,
                                              d_.W_s, d_.Ri_s);
    }
    else
    {
      update_weights_and_closest_rotations<3>(d_.Ji_m, d_.slim_energy,
                                              d_.W_m, d_.Ri_m);
      update_weights_and_closest_rotations<3>(d_.Ji_s, d_.scaf_energy,
                                              d_.W_s, d_.Ri_s);
    }
    //  cout << "update_weigths = "<<timer.getElapsedTime()<<endl;
    solve_weighted_arap(d_, uv_new);
  }
  auto whole_E =
      [&d_](Eigen::MatrixXd &uv) { return compute_energy(d_, uv); };

  igl::Timer timer;
  timer.start();
  Eigen::MatrixXi w_T;
  if (d_.m_T.cols() == d_.s_T.cols())
    igl::cat(1, d_.m_T, d_.s_T, w_T);
  else
    w_T = d_.s_T;
  double energy = igl::flip_avoiding_line_search(w_T, w_uv, V_out,
                                                 whole_E, -1) /
                  d_.mesh_measure;
  return energy;
}
}
}

IGL_INLINE void igl::scaf_precompute(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXd &V_init,
    igl::SCAFData &data,
    SLIMData::SLIM_ENERGY slim_energy,
    Eigen::VectorXi& b,
    Eigen::MatrixXd& bc,
    double soft_p)
{
  Eigen::MatrixXd CN;
  Eigen::MatrixXi FN;
  if (V_init.size() == 0)
    igl::scaf::add_new_patch(data, V, F, Eigen::RowVector2d(0, 0));
  else
    igl::scaf::add_new_patch(data, V, F, Eigen::RowVector2d(0, 0),V_init);
  data.soft_const_p = soft_p;
  igl::scaf::add_soft_constraints(data, b, bc);
  data.slim_energy = slim_energy;

  using namespace Eigen;

  auto &d_ = data;
  auto &Dx_m = data.Dx_m;
  auto &Dy_m = data.Dy_m;
  auto &Dz_m = data.Dz_m;
  auto &Dx_s = data.Dx_s;
  auto &Dy_s = data.Dy_s;
  auto &Dz_s = data.Dz_s;
  auto &Ri_m = data.Ri_m;
  auto &Ji_m = data.Ji_m;
  auto &Ri_s = data.Ri_s;
  auto &Ji_s = data.Ji_s;
  auto &W_m = data.W_m;
  auto &W_s = data.W_s;

  if (!data.has_pre_calc)
  {
    int mv_n = d_.mv_num;
    int mf_n = d_.mf_num;
    int sv_n = d_.sv_num;
    int sf_n = d_.sf_num;

    int v_n = mv_n + sv_n;
    int f_n = mf_n + sf_n;
    if (d_.dim == 2)
    {
      Eigen::MatrixXd F1, F2, F3;
      igl::local_basis(d_.m_V, d_.m_T, F1, F2, F3);
      igl::scaf::compute_surface_gradient_matrix(d_.m_V, d_.m_T, F1, F2, Dx_m,
                                                 Dy_m);

      igl::scaf::compute_scaffold_gradient_matrix(d_, Dx_s, Dy_s);
    }
    else
    {

      if (d_.m_T.cols() == 3)
      {
        igl::scaf::simplified_covariance_scatter_matrix(d_.m_V, d_.m_T,
                                                        Dx_m, Dy_m, Dz_m);
      }
      else
      {
        Eigen::SparseMatrix<double> Gm;
        igl::grad(d_.m_V, d_.m_T, Gm);

        Dx_m = Gm.block(0, 0, mf_n, mv_n);
        Dy_m = Gm.block(mf_n, 0, mf_n, mv_n);
        Dz_m = Gm.block(2 * mf_n, 0, mf_n, mv_n);
      }

      Eigen::SparseMatrix<double> Gs;
      igl::grad(d_.w_uv, d_.s_T, Gs);

      Dx_s = Gs.block(0, 0, sf_n, v_n);
      Dy_s = Gs.block(sf_n, 0, sf_n, v_n);
      Dz_s = Gs.block(2 * sf_n, 0, sf_n, v_n);
    }
    int dim = d_.dim;

    Dx_m.makeCompressed();
    Dy_m.makeCompressed();
    Dz_m.makeCompressed();
    Ri_m = MatrixXd::Zero(Dx_m.rows(), dim * dim);
    Ji_m.resize(Dx_m.rows(), dim * dim);
    W_m.resize(Dx_m.rows(), dim * dim);

    Dx_s.makeCompressed();
    Dy_s.makeCompressed();
    Dz_s.makeCompressed();
    Ri_s = MatrixXd::Zero(Dx_s.rows(), dim * dim);
    Ji_s.resize(Dx_s.rows(), dim * dim);
    W_s.resize(Dx_s.rows(), dim * dim);

    data.has_pre_calc = true;
  }
}

IGL_INLINE Eigen::MatrixXd igl::scaf_solve(SCAFData &d_, int iter_num)
{
  Eigen::VectorXi hard_cons;
  return igl::scaf_solve(d_, iter_num, hard_cons);
}

IGL_INLINE Eigen::MatrixXd igl::scaf_solve(SCAFData &d_, int iter_num, Eigen::VectorXi& cstrs)
{
  using namespace std;
  using namespace Eigen;
  double last_mesh_energy = igl::scaf::compute_energy(d_, false) / d_.mesh_measure;

  std::cout << "Initial Energy" << last_mesh_energy << std::endl;
  cout << "Initial V_num: " << d_.mv_num << " F_num: " << d_.mf_num << endl;
  for (int it = 0; it < iter_num; it++)
  {
    d_.energy = igl::scaf::compute_energy(d_, true) / d_.mesh_measure;

    igl::Timer timer;
    timer.start();
    d_.rect_frame_V = Eigen::MatrixXd();
    igl::scaf::mesh_improve(d_);

    d_.fixed_ids = cstrs;
    double new_weight = d_.mesh_measure * last_mesh_energy / (d_.sf_num * 100);
    igl::scaf::adjust_scaf_weight(d_, new_weight);

    d_.energy = igl::scaf::perform_iteration(d_);

    cout << "Iteration time = " << timer.getElapsedTime() << endl;
    double current_mesh_energy =
        igl::scaf::compute_energy(d_, false) / d_.mesh_measure;
    double mesh_energy_decrease = last_mesh_energy - current_mesh_energy;
    cout << "Energy After:"
         << d_.energy
         << "\tMesh Energy:"
         << current_mesh_energy
         << "\tEnergy Decrease"
         << mesh_energy_decrease
         << endl;
    cout << "V_num: " << d_.v_num << " F_num: " << d_.f_num << endl;
    last_mesh_energy = current_mesh_energy;
  }
  //Eigen::MatrixXd wuv3 = Eigen::MatrixXd::Zero(d_.v_num, 3);
  //wuv3.leftCols(2) = d_.w_uv;
  return d_.w_uv;
}
